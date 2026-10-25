// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_split_write.cpp
 *
 * Hadrons application, two runs: MUtilities::A2ACoarseGrid with ioSplit
 * against the same module without it (the full-grid write), one arm per
 * --io-split layout, then a read-back comparison on the coarse grid.
 *
 * Write (default), on the fine grid. Random fine vectors are sparsened under
 * the production blocking [4, 4, 4, 1], offsets 0, by every arm, each to its
 * own stem under --root: <root>/full and <root>/split_<layout>. Every arm
 * sparsens the same fine array, so every arm writes the same coarse fields.
 *
 * Read (--read), on the coarse grid: --grid is the fine one divided by the
 * blocking, with the same --mpi. --io-split here only names the stems the
 * write run produced; every stem is read the same way, on the full grid, as
 * a low-mode block of MIO::LoadCombinedA2AVecsV, and MUtilities::
 * CompareA2AVecs checks each split-written stem against the full-written
 * one. Pass is exact equality: Grid_split only copies, so the payloads are
 * identical. The files themselves are not -- the stored field norm is a
 * reduction summed in a different order on the split grid -- which is why
 * the check is on the read-back fields rather than the bytes. The read checks
 * each stored norm against the field, so the norm record is exercised too.
 *
 * Pick --size so it is not a multiple of the partition count, so the last
 * batch has idle partitions. Clear --root between write runs.
 *
 * Usage:
 *   srun ... ./Test_split_write --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --root /path/split_write --io-split 2.2.2.1,2.2.2.4
 *   srun ... ./Test_split_write --grid 16.16.16.128 --mpi 4.4.4.8 \
 *       --root /path/split_write --io-split 2.2.2.1,2.2.2.4 --read
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

using namespace Grid;
using namespace Hadrons;

static const std::vector<int> BLOCK_SIZE = {4, 4, 4, 1};
static const std::vector<int> OFFSETS    = {0, 0, 0, 0};

static void addSparsen(Application &application, const std::string &name,
                       const std::string &fine, const std::string &output,
                       const std::string &ioSplit)
{
    MUtilities::A2ACoarseGridPar par;

    par.fine       = fine;
    par.modeOffset = 0;
    par.modeSize   = 0;
    par.scale      = 1.0;
    par.blockSize  = BLOCK_SIZE;
    par.offsets   = OFFSETS;
    par.output    = output;
    par.ioSplit   = ioSplit;

    application.createModule<MUtilities::A2ACoarseGrid>(name, par);
}

static void addLoad(Application &application, const std::string &name,
                    const std::string &stem, const unsigned int size)
{
    MIO::LoadCombinedA2AVecsVPar par;

    par.lowFilestem = stem;
    par.nLow        = size;
    par.highStem    = "";
    par.highSize    = 0;
    par.nHit        = 0;

    application.createModule<MIO::LoadCombinedA2AVecsV>(name, par);
}

static void addCompare(Application &application, const std::string &name,
                       const std::string &left, const std::string &right)
{
    MUtilities::CompareA2AVecsPar par;

    par.left  = left;
    par.right = right;

    application.createModule<MUtilities::CompareA2AVecs>(name, par);
}

// "2.2.2.1" -> "2 2 2 1", the form strToVec reads.
static std::string layoutToPar(std::string layout)
{
    for (auto &c : layout)
        if (c == '.') c = ' ';

    return layout;
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
    int          traj    = 0;
    unsigned int size    = 84;
    bool         read    = GridCmdOptionExists(argv, argv + argc, "--read");
    std::string  root, splitList;

    if (GridCmdOptionExists(argv, argv + argc, "--traj"))
        traj      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--traj"));
    if (GridCmdOptionExists(argv, argv + argc, "--size"))
        size      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--size"));
    if (GridCmdOptionExists(argv, argv + argc, "--root"))
        root      = GridCmdOptionPayload(argv, argv + argc, "--root");
    if (GridCmdOptionExists(argv, argv + argc, "--io-split"))
        splitList = GridCmdOptionPayload(argv, argv + argc, "--io-split");

    std::vector<std::string> layouts = strToVec<std::string>(cliListToPar(splitList));

    if (root.empty())
    {
        HADRONS_ERROR(Argument, "--root is required");
    }
    if (layouts.empty())
    {
        HADRONS_ERROR(Argument, "--io-split needs at least one layout");
    }

    // ------------------------------------------------------------------
    // Global parameters. Naive scheduling keeps the arms in creation order.
    // ------------------------------------------------------------------
    Application application;
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start       = traj;
    globalPar.trajCounter.end         = traj + 1;
    globalPar.trajCounter.step        = 1;
    globalPar.runId                   = read ? "split_write_read" : "split_write";
    globalPar.scheduler.schedulerType = "naive";
    globalPar.genetic.maxGen          = 1000;
    globalPar.genetic.maxCstGen       = 200;
    globalPar.genetic.popSize         = 20;
    globalPar.genetic.mutationRate    = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // The arms.
    // ------------------------------------------------------------------
    if (read)
    {
        addLoad(application, "v_written_full", root + "/full", size);
        for (auto &l : layouts)
        {
            addLoad(application, "v_written_split_" + l, root + "/split_" + l, size);
            addCompare(application, "cmp_" + l, "v_written_full",
                       "v_written_split_" + l);
        }
    }
    else
    {
        MUtilities::RandomVectorsPar rv;
        rv.size      = size;
        rv.Ls        = 1;
        rv.output    = "";
        rv.multiFile = false;
        application.createModule<MUtilities::RandomFermions>("fine", rv);

        addSparsen(application, "sp_full", "fine", root + "/full", "");
        for (auto &l : layouts)
        {
            addSparsen(application, "sp_split_" + l, "fine", root + "/split_" + l,
                       layoutToPar(l));
        }
    }

    application.run();

    Grid_finalize();

    return EXIT_SUCCESS;
}
