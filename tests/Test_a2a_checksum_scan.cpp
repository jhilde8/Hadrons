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
 * Usage:
 *   srun ... ./Test_a2a_checksum_scan --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --io-split 2.2.2.1 --traj 1540 --first 1024 --size 512 \
 *       --stems /path/vw/sloppy/s8_v,/path/vw/exact/s8_exact
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
    template <class vobj>
    int verify(Lattice<vobj> &field, RealD &rdiff)
    {
        typedef typename vobj::scalar_object sobj;

        std::string       record(ILDG_BINARY_DATA);
        scidacChecksum    checksum;
        FieldNormMetaData normMeta;
        uint32_t          nersc, suma, sumb;
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
    std::string  stemList, ioSplit;

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

    std::vector<std::string> stems = strToVec<std::string>(cliListToPar(stemList));

    if (stems.empty())
    {
        HADRONS_ERROR(Argument, "--stems needs at least one stem");
    }
    if (size == 0)
    {
        HADRONS_ERROR(Argument, "--size is required");
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

    unsigned int nPart     = grid->ProcessorCount()/ioGrid->ProcessorCount();
    unsigned int nBadTotal = 0;
    Field        field(ioGrid);

    LOG(Message) << "Scanning " << stems.size() << " stem(s), files " << first
                 << " - " << first + size - 1 << ", on " << nPart
                 << " partition(s)" << std::endl;

    // ------------------------------------------------------------------
    // Scan.
    // ------------------------------------------------------------------
    for (auto &stem: stems)
    {
        std::vector<RealD> status(size, 0.), rdiff(size, 0.);

        for (unsigned int b = 0; b < size; b += nPart)
        {
            unsigned int i = b + part;

            if (i < size)
            {
                std::string filename = elementFilename(stem, traj, first + i);
                int         s        = Missing;
                RealD       r        = 0.;

                if (std::ifstream(filename).good())
                {
                    ScanReader reader;

                    reader.open(filename);
                    s = reader.verify(field, r);
                    reader.close();
                }
                if (ioGrid->IsBoss())
                {
                    status[i] = s;
                    rdiff[i]  = r;
                }
            }
            LOG(Message) << stem << ": files " << first + b << " - "
                         << first + std::min(b + nPart, size) - 1 << " read" << std::endl;
        }
        grid->GlobalSumVector(status.data(), size);
        grid->GlobalSumVector(rdiff.data(), size);

        unsigned int nBad = 0;

        for (unsigned int i = 0; i < size; ++i)
        {
            int s = static_cast<int>(status[i]);

            if (s != Good)
            {
                LOG(Message) << "BAD " << elementFilename(stem, traj, first + i)
                             << ":" << statusString(s) << " (norm rdiff "
                             << rdiff[i] << ")" << std::endl;
                nBad++;
            }
        }
        LOG(Message) << stem << ": " << nBad << " bad of " << size << std::endl;
        nBadTotal += nBad;
    }
    LOG(Message) << "Total: " << nBadTotal << " bad" << std::endl;

    Grid_finalize();

    return (nBadTotal == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
