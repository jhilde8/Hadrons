/*
 * Test_a2a_io_width.cpp
 *
 * A2AVectorsIo write/read bandwidth versus SciDAC record width, to settle the
 * on-disk format for the A2A vectors once binning is retired.
 *
 * One mode per run, all four record types inside it, so a single job yields the
 * comparison under one set of filesystem conditions:
 *
 *   fermion   FIMPL::FermionField                    192 B/site
 *   bin1      Lattice<iVector<SiteSpinor,1>>         192 B/site
 *   prop      FIMPL::PropagatorField                2304 B/site
 *   bin12     Lattice<iVector<SiteSpinor,12>>       2304 B/site
 *
 * bin12 is byte-for-byte the current production format, prop is the same width
 * through a Grid-native type, and fermion is the width-1 candidate. --nvec
 * counts fermion-equivalent vectors and is divided by 12 for the wide types, so
 * every type moves the SAME total bytes: the comparison is record width, not
 * volume. Each type gets its own stem (<stem>_fermion, <stem>_bin1, ...) and
 * its own subdirectory, so nothing collides.
 *
 * The two pairs exist because byte width is not the only variable. prop and
 * bin12 are the same size per site but a different tensor nest, so BinaryIO's
 * scalar-object gather, checksum and byte-order passes walk them differently;
 * fermion and bin1 are the same pair at width 1. A gap inside a pair is layout
 * or filesystem drift, never width, which is what makes them the control on a
 * gap between pairs. Measured 2026-10-08, write at 4 nodes / 32 ranks on the
 * 64I volume (77.3 GB per type): fermion 3276 MB/s, bin12 2331, prop 2025 --
 * so the narrow record was faster, and two types of equal width differed by
 * 15%.
 *
 * That 15% was taken with one record for each wide type, which is also the one
 * configuration where the first-record cost cannot be separated out: the four
 * BinaryIO host buffers are allocated and first-touched on record 0 and then
 * come back from Grid's allocator cache, so with nrec = 1 that cost sits
 * entirely inside the only sample. The summary therefore reports record 0
 * separately from the mean of the rest. Give the wide types at least three or
 * four records (--nvec 36 or 48) before reading anything into a gap.
 *
 * Everything goes through A2AVectorsIo::writeElement / readElement, which is
 * one file per record at <stem>_<type>.<traj>/elem<N>.bin -- the layout the
 * post-binning format would use, and already what MIO::LoadBinnedA2AVecs1
 * (fermion) and MIO::LoadBinnedA2AVecs12 (bin12) read, so a module-level
 * cross-check needs no new registration.
 *
 * Write and read are separate runs because in production they happen at
 * different node counts: the CG writes, the contractions read, and with the A2A
 * loop and sparsening possibly moving to their own smaller job the two ends
 * decouple further. There is deliberately no combined mode -- reading back what
 * this process just wrote measures client cache, not the filesystem.
 *
 * Peak memory is one record field (two with --verify), not the whole set, so
 * either mode runs at any node count.
 *
 * --verify, in read mode, regenerates each source from the seed and compares.
 * It assumes SeedFixedIntegers is independent of the MPI decomposition, which
 * is what makes it usable when the read job has a different --mpi from the
 * write job; it needs the same --seed and --nvec as the write run.
 *
 * What to read off the output: the summary table, and Grid's own
 *   IOobject: aggregate N extent(s)/rank, first X MB (target 4 MB)
 * lines. If X stays near 4 MB for fermion, the aggregator is simply gathering
 * more extents to reach the same byte target and the narrow record costs
 * nothing; if it collapses, the wide record is buying something real. Sweeping
 * node count at fixed --nvec maps bandwidth against bytes/rank, which is the
 * number that predicts the sparsened-field case without building one.
 *
 * The three types are written as consecutive blocks, so filesystem drift lands
 * on whichever block is running. The per-record min/max columns expose that;
 * if the gaps between types are of that order, run it again before believing
 * them.
 *
 * Stripe the output directory first -- a 6-77 GB file sits entirely inside the
 * OLCF PFL default's -c 1 component:
 *   lfs setstripe -c 16 <dir>
 *
 * Usage:
 *   srun -n 512 ./Test_a2a_io_width --grid 64.64.64.128 --mpi 4.4.4.8 \
 *        --mode write --nvec 24 --stem $SCRATCH/iotest/v
 *   srun -n 4096 ./Test_a2a_io_width --grid 64.64.64.128 --mpi 8.8.8.8 \
 *        --mode read --nvec 24 --stem $SCRATCH/iotest/v --verify
 *   mpirun -n 1 ./Test_a2a_io_width --grid 8.8.8.16 --mpi 1.1.1.1 \
 *        --mode write --nvec 12 --record bin12 --stem ./iotest/v
 *
 * A 32.32.32.32 volume on 4.4.2.4 reproduces the per-rank geometry of the 64I
 * volume on 8.8.8.8 exactly -- 8192 sites/rank, an 8-site contiguous extent,
 * 1024 extents/rank -- at 1/32 of the ranks, so the extent-aggregation question
 * can be answered at 16 nodes. Absolute bandwidth, per-file open cost and the
 * ROMIO aggregator topology do not carry over from that.
 */

#include <iomanip>
#include <sstream>
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AVectors.hpp>

using namespace Grid;
using namespace Hadrons;

namespace
{

struct IoResult
{
    std::string  label;
    unsigned int bytesPerSite{0};
    double       bytesPerRec{0.};
    double       bytesPerRankPerRec{0.};
    unsigned int nRec{0};
    double       total{0.};
    double       min{0.};
    double       max{0.};
    // Record 0 carries the allocation and first-touch of BinaryIO's host
    // buffers, which later records get back from Grid's allocator cache.
    double       first{0.};
    // Negative when no verification was requested or possible.
    double       maxDiff{-1.};
};

// One record type end to end. The RNG is seeded at the top of each pass so a
// read run can regenerate exactly the fields a write run stored, without
// either pass holding more than one record resident.
template <typename Field>
IoResult runIo(const std::string label, const std::string stem,
               const unsigned int nRec, const int traj, const bool doWrite,
               const bool verify, GridCartesian *grid,
               const std::vector<int> &seed)
{
    typedef typename Field::vector_object vobj;
    typedef typename vobj::scalar_object  sobj;

    IoResult        r;
    GridParallelRNG rng(grid);
    Field           f(grid);
    double          tMin = std::numeric_limits<double>::max(), tMax = 0.;

    r.label              = label;
    r.bytesPerSite       = sizeof(sobj);
    r.bytesPerRec        = (double)sizeof(sobj)*(double)grid->gSites();
    r.bytesPerRankPerRec = r.bytesPerRec/(double)grid->ProcessorCount();
    r.nRec               = nRec;

    LOG(Message) << "=== " << label << ": " << r.bytesPerSite << " B/site, "
                 << r.bytesPerRec/1.0e9 << " GB/record, " << nRec
                 << " record(s), " << r.bytesPerRankPerRec/1.0e6
                 << " MB/rank/record ===" << std::endl;

    rng.SeedFixedIntegers(seed);
    if (doWrite)
    {
        for (unsigned int i = 0; i < nRec; ++i)
        {
            double t;

            random(rng, f);
            t = -usecond();
            A2AVectorsIo::writeElement(stem, f, i, traj);
            t += usecond();
            t /= 1.0e6;
            r.total += t;
            tMin     = std::min(tMin, t);
            tMax     = std::max(tMax, t);
            if (i == 0) { r.first = t; }
        }
    }
    else
    {
        Field ref(grid);

        for (unsigned int i = 0; i < nRec; ++i)
        {
            double t;

            t = -usecond();
            A2AVectorsIo::readElement(stem, f, i, traj);
            t += usecond();
            t /= 1.0e6;
            r.total += t;
            tMin     = std::min(tMin, t);
            tMax     = std::max(tMax, t);
            if (i == 0) { r.first = t; }
            if (verify)
            {
                random(rng, ref);
                ref       -= f;
                r.maxDiff  = std::max(std::max(r.maxDiff, 0.), norm2(ref));
            }
        }
    }
    r.min = (nRec > 0) ? tMin : 0.;
    r.max = tMax;

    LOG(Message) << label << ": " << r.total << " s" << std::endl;

    return r;
}

void summary(const std::string mode, const std::vector<IoResult> &res)
{
    LOG(Message) << std::endl;
    LOG(Message) << "================ " << mode << " summary ================"
                 << std::endl;
    LOG(Message) << std::left
                 << std::setw(10) << "record"
                 << std::right
                 << std::setw(9)  << "B/site"
                 << std::setw(11) << "GB/rec"
                 << std::setw(13) << "MB/rank/rec"
                 << std::setw(7)  << "nrec"
                 << std::setw(11) << "total GB"
                 << std::setw(11) << "total s"
                 << std::setw(11) << "MB/s"
                 << std::setw(11) << "s/rec"
                 << std::endl;
    for (auto &r: res)
    {
        const double totalBytes = r.bytesPerRec*r.nRec;

        LOG(Message) << std::left
                     << std::setw(10) << r.label
                     << std::right << std::fixed
                     << std::setw(9)  << r.bytesPerSite
                     << std::setprecision(3)
                     << std::setw(11) << r.bytesPerRec/1.0e9
                     << std::setw(13) << r.bytesPerRankPerRec/1.0e6
                     << std::setw(7)  << r.nRec
                     << std::setw(11) << totalBytes/1.0e9
                     << std::setw(11) << r.total
                     << std::setprecision(1)
                     << std::setw(11) << ((r.total > 0.) ? (totalBytes/1.0e6)/r.total : 0.)
                     << std::setprecision(3)
                     << std::setw(11) << ((r.nRec > 0) ? r.total/r.nRec : 0.)
                     << std::endl;
    }
    LOG(Message) << "per record, record 0 separated (it carries the buffer "
                 << "allocation and first touch):" << std::endl;
    for (auto &r: res)
    {
        std::stringstream rest;

        rest << std::fixed << std::setprecision(3);
        if (r.nRec > 1)
        {
            rest << (r.total - r.first)/(r.nRec - 1) << " s";
        }
        else
        {
            rest << "n/a (single record, raise --nvec)";
        }
        LOG(Message) << "  " << std::left << std::setw(10) << r.label
                     << std::right << std::fixed << std::setprecision(3)
                     << " first " << r.first << " s, rest mean " << rest.str()
                     << ", min " << r.min << " s, max " << r.max << " s"
                     << std::endl;
    }
    for (auto &r: res)
    {
        if (r.maxDiff >= 0.)
        {
            LOG(Message) << "verify " << r.label << ": max norm2 difference "
                         << r.maxDiff << std::endl;
            if (r.maxDiff != 0.)
            {
                LOG(Error) << "round trip is not exact for record type '"
                           << r.label << "'" << std::endl;
            }
        }
    }
}

}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    std::string      mode   = "";
    std::string      record = "all";
    std::string      stem   = "iotest/v";
    unsigned int     nVec   = 24;
    int              traj   = 0;
    bool             verify = false;
    std::vector<int> seed   = {1, 2, 3, 4};

    if (GridCmdOptionExists(argv, argv + argc, "--mode"))
        mode   = GridCmdOptionPayload(argv, argv + argc, "--mode");
    if (GridCmdOptionExists(argv, argv + argc, "--record"))
        record = GridCmdOptionPayload(argv, argv + argc, "--record");
    if (GridCmdOptionExists(argv, argv + argc, "--stem"))
        stem   = GridCmdOptionPayload(argv, argv + argc, "--stem");
    if (GridCmdOptionExists(argv, argv + argc, "--nvec"))
        nVec   = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--nvec"));
    if (GridCmdOptionExists(argv, argv + argc, "--traj"))
        traj   = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--traj"));
    if (GridCmdOptionExists(argv, argv + argc, "--verify"))
        verify = true;
    if (GridCmdOptionExists(argv, argv + argc, "--seed"))
        GridCmdOptionIntVector(GridCmdOptionPayload(argv, argv + argc, "--seed"), seed);

    const bool doWrite = (mode == "write");
    const bool doRead  = (mode == "read");
    const bool doFerm  = (record == "all") || (record == "fermion");
    const bool doBin1  = (record == "all") || (record == "bin1");
    const bool doProp  = (record == "all") || (record == "prop");
    const bool doBin12 = (record == "all") || (record == "bin12");

    if (!doWrite && !doRead)
    {
        LOG(Error) << "--mode must be given explicitly as write or read"
                   << (mode.empty() ? "" : (" (got '" + mode + "')"))
                   << "; there is no combined mode, see the file header"
                   << std::endl;
        Grid_finalize();
        return EXIT_FAILURE;
    }
    if (!doFerm && !doBin1 && !doProp && !doBin12)
    {
        LOG(Error) << "--record must be all, fermion, bin1, prop or bin12 "
                   << "(got '" << record << "')" << std::endl;
        Grid_finalize();
        return EXIT_FAILURE;
    }
    if (nVec == 0)
    {
        LOG(Error) << "--nvec must be positive" << std::endl;
        Grid_finalize();
        return EXIT_FAILURE;
    }
    if ((doProp || doBin12) && (nVec % 12 != 0))
    {
        LOG(Error) << "--nvec must be a multiple of 12 so that every record "
                   << "type moves the same total bytes (got " << nVec << ")"
                   << std::endl;
        Grid_finalize();
        return EXIT_FAILURE;
    }
    if (verify && doWrite)
    {
        LOG(Warning) << "--verify only applies to --mode read; ignoring it"
                     << std::endl;
        verify = false;
    }

    GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi());

    LOG(Message) << "A2AVectorsIo record width test" << std::endl;
    LOG(Message) << "mode '" << mode << "', record '" << record << "', nvec "
                 << nVec << ", traj " << traj << ", verify " << verify
                 << std::endl;
    LOG(Message) << "stem '" << stem << "' -> " << stem << "_<type>." << traj
                 << "/elem<N>.bin" << std::endl;

    std::vector<IoResult> res;

    if (doFerm)
    {
        res.push_back(runIo<FIMPL::FermionField>(
            "fermion", stem + "_fermion", nVec, traj, doWrite, verify, grid,
            seed));
    }
    if (doBin1)
    {
        typedef Lattice<iVector<FIMPL::SiteSpinor, 1>> Bin1Field;

        res.push_back(runIo<Bin1Field>(
            "bin1", stem + "_bin1", nVec, traj, doWrite, verify, grid, seed));
    }
    if (doProp)
    {
        res.push_back(runIo<FIMPL::PropagatorField>(
            "prop", stem + "_prop", nVec/12, traj, doWrite, verify, grid,
            seed));
    }
    if (doBin12)
    {
        typedef Lattice<iVector<FIMPL::SiteSpinor, 12>> Bin12Field;

        res.push_back(runIo<Bin12Field>(
            "bin12", stem + "_bin12", nVec/12, traj, doWrite, verify, grid,
            seed));
    }

    summary(doWrite ? "WRITE" : "READ", res);

    Grid_finalize();

    return EXIT_SUCCESS;
}
