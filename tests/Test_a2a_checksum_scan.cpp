/*
 * Test_a2a_checksum_scan.cpp
 *
 * Checks every A2A vector file of one or more stems against its stored SciDAC
 * checksums and field norm, and lists every bad file instead of aborting on
 * the first one as the loaders do. A file is bad if either checksum differs,
 * if its norm differs from the stored one by rdiff >= 1e-5 (the loaders'
 * GRID_FIELD_NORM_CHECK), or if it is missing or has no binary record.
 *
 * A read error (EIO) still aborts: it is a GRID_ASSERT inside BinaryIO, and
 * Lustre returns it on a transient server timeout as well as on a damaged
 * file. Rerun from the stem that failed; a file that fails every time is
 * damaged.
 *
 * The read is ScanReader::verify, which is GridLimeReader::
 * readLimeLatticeBinaryObject with the checksum and norm asserts turned into
 * a status. With --io-split each partition reads a different file of a batch,
 * as in A2AVectorsIo::readBlock, but nothing is unsplit: each partition's
 * boss records its verdict and one GlobalSumVector per stem gathers them.
 *
 * --grid is the grid the files were written on (the coarse one for sparsened
 * vectors). Files are elem<first> .. elem<first+size-1> of <stem>.<traj>/
 * (--first defaults to 0), the same range for every stem, so the low modes and
 * the hits are separate runs. The exit code is nonzero if any file is bad.
 *
 * Replica comparison (--compare a,b): on a machine that can return wrong
 * answers every vector job runs twice, each replica writing its own tree
 * (PyXmlGen config.vw_base, vw.<tag>). Every stem then holds the placeholder
 * @replica@, replaced by each tag in turn, and every replica's files are
 * scanned as above. A wrong solve writes a file that agrees with its own
 * checksum, so the scan alone cannot see it; but its data differs from the
 * other replica's, and so does the SciDAC checksum pair computed from it. An
 * element whose files are good in every replica is compared on that pair
 * against the first replica: equal pairs are bitwise agreement up to a
 * 2^-64 collision. Any difference is a MISMATCH, also a nonzero exit; an
 * element bad in some replica cannot be compared and is reported as such.
 *
 * Usage:
 *   srun ... ./Test_a2a_checksum_scan --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --io-split 2.2.2.1 --traj 1540 --first 1024 --size 512 \
 *       --stems /path/vw/sloppy/s8_v,/path/vw/exact/s8_exact
 *   mpiexec ... ./Test_a2a_checksum_scan ... --traj 1200 --size 1536 \
 *       --stems /path/vw.@replica@/sloppy/s3_v --compare a,b
 */

#include <Hadrons/Global.hpp>
#include <Hadrons/A2AVectors.hpp>
#include "TestMFShells.hpp"

using namespace Grid;
using namespace Hadrons;

typedef FIMPL::FermionField Field;

enum ScanStatus
{
    Good        = 0,
    BadChecksum = 1,
    BadNorm     = 2,
    Missing     = 4
};

class ScanReader: public ScidacReader
{
public:
    // suma / sumb are the checksums computed from the data read, which equal
    // the stored ones when the status is Good.
    template <class vobj>
    int verify(Lattice<vobj> &field, RealD &rdiff, uint32_t &suma, uint32_t &sumb)
    {
        typedef typename vobj::scalar_object sobj;

        std::string       record(ILDG_BINARY_DATA);
        scidacChecksum    checksum;
        FieldNormMetaData normMeta;
        uint32_t          nersc;
        int               status = Good;

        while (limeReaderNextRecord(LimeR) == LIME_SUCCESS)
        {
            if (!strncmp(limeReaderType(LimeR), record.c_str(), record.size()))
            {
                GRID_ASSERT(sizeof(sobj)*field.Grid()->_gsites == limeReaderBytes(LimeR));

                uint64_t                       offset = ftello(File);
                BinarySimpleMunger<sobj, sobj> munge;

                BinaryIO::readLatticeObject<vobj, sobj>(field, filename, munge, offset,
                                                        getFormatString<vobj>(),
                                                        nersc, suma, sumb,
                                                        BINARYIO_LEXICOGRAPHIC|BINARYIO_AGGREGATE);
                readScidacChecksum(checksum, normMeta);
                if (normMeta.norm2 != 0.0)
                {
                    rdiff = GRID_FIELD_NORM_CALC(normMeta, norm2(field));
                    if (rdiff >= 1.0e-5)
                    {
                        status |= BadNorm;
                    }
                }
                if (scidacChecksumVerify(checksum, suma, sumb) != 1)
                {
                    status |= BadChecksum;
                }

                return status;
            }
        }

        return Missing;
    }
};

// A2AVectorsIo::elementFilename, which is private.
static std::string elementFilename(const std::string &stem, const int traj,
                                   const unsigned int index)
{
    std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

    return stem + t + "/elem" + std::to_string(index) + ".bin";
}

static const std::string replicaKey = "@replica@";

static std::string replicaStem(const std::string &stem, const std::string &tag)
{
    std::string s   = stem;
    size_t      pos = s.find(replicaKey);

    if (pos != std::string::npos)
    {
        s.replace(pos, replicaKey.size(), tag);
    }

    return s;
}

// One replica's scan of one stem, indexed by file. Each entry is set by the
// boss of the partition that read the file and zero elsewhere, so one
// GlobalSumVector per array gathers them; checksums are exact in a RealD.
struct StemScan
{
    std::vector<RealD> status, rdiff, suma, sumb;

    explicit StemScan(const unsigned int size)
    : status(size, 0.), rdiff(size, 0.), suma(size, 0.), sumb(size, 0.)
    {}

    void gather(GridBase *grid)
    {
        unsigned int size = status.size();

        grid->GlobalSumVector(status.data(), size);
        grid->GlobalSumVector(rdiff.data(), size);
        grid->GlobalSumVector(suma.data(), size);
        grid->GlobalSumVector(sumb.data(), size);
    }
};

static std::string checksumString(const StemScan &scan, const unsigned int i)
{
    std::ostringstream s;

    s << std::hex << static_cast<uint32_t>(scan.suma[i]) << " "
      << static_cast<uint32_t>(scan.sumb[i]);

    return s.str();
}

static std::string statusString(const int status)
{
    std::string s;

    if (status & BadChecksum) s += " checksum";
    if (status & BadNorm)     s += " norm";
    if (status & Missing)     s += " missing";

    return s;
}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    // ------------------------------------------------------------------
    // Parse CLI arguments.
    // ------------------------------------------------------------------
    int          traj = 0;
    unsigned int first = 0, size = 0;
    std::string  stemList, ioSplit, compareList;

    if (GridCmdOptionExists(argv, argv + argc, "--traj"))
        traj     = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--traj"));
    if (GridCmdOptionExists(argv, argv + argc, "--first"))
        first    = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--first"));
    if (GridCmdOptionExists(argv, argv + argc, "--size"))
        size     = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--size"));
    if (GridCmdOptionExists(argv, argv + argc, "--stems"))
        stemList = GridCmdOptionPayload(argv, argv + argc, "--stems");
    if (GridCmdOptionExists(argv, argv + argc, "--io-split"))
        ioSplit  = GridCmdOptionPayload(argv, argv + argc, "--io-split");
    if (GridCmdOptionExists(argv, argv + argc, "--compare"))
        compareList = GridCmdOptionPayload(argv, argv + argc, "--compare");

    std::vector<std::string> stems    = strToVec<std::string>(cliListToPar(stemList));
    std::vector<std::string> replicas = strToVec<std::string>(cliListToPar(compareList));

    if (stems.empty())
    {
        HADRONS_ERROR(Argument, "--stems needs at least one stem");
    }
    if (size == 0)
    {
        HADRONS_ERROR(Argument, "--size is required");
    }
    if (!compareList.empty() && replicas.size() < 2)
    {
        HADRONS_ERROR(Argument, "--compare needs at least two replica tags");
    }
    for (auto &stem: stems)
    {
        if ((stem.find(replicaKey) != std::string::npos) != !replicas.empty())
        {
            HADRONS_ERROR(Argument, "stem '" + stem + "': " + replicaKey
                          + " is required with --compare and not allowed without");
        }
    }
    // A plain scan is a comparison of one unnamed replica with nothing
    if (replicas.empty())
    {
        replicas.push_back("");
    }

    // ------------------------------------------------------------------
    // Grids. Without --io-split the full grid reads every file in turn.
    // ------------------------------------------------------------------
    GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()), GridDefaultMpi());
    std::unique_ptr<GridCartesian> splitGrid;
    GridBase                       *ioGrid = grid;
    unsigned int                   part    = 0;

    if (!ioSplit.empty())
    {
        for (auto &c : ioSplit)
            if (c == '.') c = ' ';
        splitGrid.reset(new GridCartesian(grid->FullDimensions(), grid->_simd_layout,
                                          strToVec<int>(ioSplit), *grid));
        ioGrid = splitGrid.get();
        part   = GridSplitVectorIndex(grid, splitGrid.get());
    }

    unsigned int nPart          = grid->ProcessorCount()/ioGrid->ProcessorCount();
    unsigned int nRep           = replicas.size();
    unsigned int nBadTotal      = 0;
    unsigned int nMismatchTotal = 0;
    Field        field(ioGrid);

    LOG(Message) << "Scanning " << stems.size() << " stem(s), files " << first
                 << " - " << first + size - 1 << ", on " << nPart
                 << " partition(s)"
                 << ((nRep > 1) ? ", replicas " + compareList + " compared" : "")
                 << std::endl;

    // ------------------------------------------------------------------
    // Scan every replica of every stem.
    // ------------------------------------------------------------------
    for (auto &stem: stems)
    {
        std::vector<StemScan> scans(nRep, StemScan(size));

        for (unsigned int r = 0; r < nRep; ++r)
        {
            std::string repStem = replicaStem(stem, replicas[r]);
            StemScan    &scan   = scans[r];

            for (unsigned int b = 0; b < size; b += nPart)
            {
                unsigned int i = b + part;

                if (i < size)
                {
                    std::string filename = elementFilename(repStem, traj, first + i);
                    int         s        = Missing;
                    RealD       rd       = 0.;
                    uint32_t    suma     = 0, sumb = 0;

                    if (std::ifstream(filename).good())
                    {
                        ScanReader reader;

                        reader.open(filename);
                        s = reader.verify(field, rd, suma, sumb);
                        reader.close();
                    }
                    if (ioGrid->IsBoss())
                    {
                        scan.status[i] = s;
                        scan.rdiff[i]  = rd;
                        scan.suma[i]   = suma;
                        scan.sumb[i]   = sumb;
                    }
                }
                LOG(Message) << repStem << ": files " << first + b << " - "
                             << first + std::min(b + nPart, size) - 1 << " read" << std::endl;
            }
            scan.gather(grid);

            unsigned int nBad = 0;

            for (unsigned int i = 0; i < size; ++i)
            {
                int s = static_cast<int>(scan.status[i]);

                if (s != Good)
                {
                    LOG(Message) << "BAD " << elementFilename(repStem, traj, first + i)
                                 << ":" << statusString(s) << " (norm rdiff "
                                 << scan.rdiff[i] << ")" << std::endl;
                    nBad++;
                }
            }
            LOG(Message) << repStem << ": " << nBad << " bad of " << size << std::endl;
            nBadTotal += nBad;
        }

        // --------------------------------------------------------------
        // Compare the replicas' checksums, element by element.
        // --------------------------------------------------------------
        if (nRep > 1)
        {
            unsigned int nMismatch = 0, nUncompared = 0;

            for (unsigned int i = 0; i < size; ++i)
            {
                bool allGood = true, agree = true;

                for (unsigned int r = 0; r < nRep; ++r)
                {
                    allGood = allGood && (static_cast<int>(scans[r].status[i]) == Good);
                    agree   = agree && (scans[r].suma[i] == scans[0].suma[i])
                                    && (scans[r].sumb[i] == scans[0].sumb[i]);
                }
                if (!allGood)
                {
                    nUncompared++;
                }
                else if (!agree)
                {
                    std::ostringstream line;

                    for (unsigned int r = 0; r < nRep; ++r)
                    {
                        line << " " << replicas[r] << " " << checksumString(scans[r], i);
                    }
                    LOG(Message) << "MISMATCH " << stem << " elem " << first + i
                                 << ":" << line.str() << std::endl;
                    nMismatch++;
                }
            }
            LOG(Message) << stem << ": " << nMismatch << " mismatch(es), "
                         << nUncompared << " not compared (bad in a replica), of "
                         << size << std::endl;
            nMismatchTotal += nMismatch;
        }
    }
    LOG(Message) << "Total: " << nBadTotal << " bad"
                 << ((nRep > 1) ? ", " + std::to_string(nMismatchTotal) + " mismatch(es)" : "")
                 << std::endl;

    Grid_finalize();

    return (nBadTotal == 0 && nMismatchTotal == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
