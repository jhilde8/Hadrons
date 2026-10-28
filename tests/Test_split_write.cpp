// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_split_write.cpp
 *
 * Hadrons application, two runs: a writer module with ioSplit against the
 * same module without it (the full-grid write), one arm per --io-split
 * layout, then a read-back comparison. Two writers:
 *
 *   default  MUtilities::A2ACoarseGrid. Random fine vectors are sparsened
 *            under the production blocking [4, 4, 4, 1], offsets 0, to
 *            <root>/full and <root>/split_<layout>.
 *   --low    MUtilities::A2ALowModesCoarse200 on the first --size coarse
 *            vectors of the 64I light eigenpack (EPACK_STEM, 200 fine), with
 *            the production light action on a unit gauge field: unphysical
 *            low modes, but the same ones in every arm. Writes
 *            <root>/full_v, <root>/full_w and <root>/split_<layout>_{v,w}.
 *            Each batch is one writeBlock call from mode il0, so every batch
 *            after the first exercises a nonzero firstIndex.
 *
 * Every arm writes the same fields from the same input.
 *
 * Read (--read, with --low again for the low modes). The sparsened fields
 * are read on the coarse grid: --grid is the fine one divided by the
 * blocking, with the same --mpi. The low modes are read on the same grid as
 * they were written. --io-split here only names the stems the write run
 * produced; every stem is read the same way, on the full grid, as a low-mode
 * block of MIO::LoadCombinedA2AVecsV, and MUtilities::CompareA2AVecs checks
 * each split-written stem against the full-written one. Pass is exact
 * equality: Grid_split only copies, so the payloads are identical. The files
 * themselves are not -- the stored field norm is a reduction summed in a
 * different order on the split grid -- which is why the check is on the
 * read-back fields rather than the bytes. The read checks each stored norm
 * against the field, so the norm record is exercised too.
 *
 * Pick --size so it is not a multiple of the partition count but larger than
 * it, so there are several batches and the last has idle partitions. A low
 * mode is a full-volume field, 6.4 GB a file at 64^3 x 128: size the run
 * against the space under --root. Clear --root between write runs.
 *
 * Usage:
 *   srun ... ./Test_split_write --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --root /path/split_write --io-split 2.2.2.1,2.2.2.4
 *   srun ... ./Test_split_write --grid 16.16.16.128 --mpi 4.4.4.8 \
 *       --root /path/split_write --io-split 2.2.2.1,2.2.2.4 --read
 *   srun ... ./Test_split_write --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --root /path/split_low --io-split 2.1.2.2,2.2.2.4 --size 70 --low
 *   srun ... ./Test_split_write --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --root /path/split_low --io-split 2.1.2.2,2.2.2.4 --size 70 --low --read
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

using namespace Grid;
using namespace Hadrons;

static const std::vector<int> BLOCK_SIZE = {4, 4, 4, 1};
static const std::vector<int> OFFSETS    = {0, 0, 0, 0};

// PyXmlGen/config.py CG_EIGENPACK and the light CG_ACTIONS entry.
static const std::string EPACK_STEM =
    "/lustre/orion/phy157/proj-shared/phy157_dwf/jhilde/main_64I/converted/vec";

static void addLowModeInputs(Application &application, const unsigned int size)
{
    MIO::LoadCoarseEigenPackPar ep;
    MAction::MobiusDWFPar       act;

    application.createModule<MGauge::Unit>("gauge");

    ep.filestem      = EPACK_STEM;
    ep.multiFile     = true;
    ep.sizeFine      = 200;
    ep.sizeCoarse    = size;
    ep.redBlack      = true;
    ep.Ls            = 12;
    ep.blockSize     = "4 4 4 4 12";
    ep.orthogonalise = false;
    ep.gaugeXform    = "";
    application.createModule<MIO::LoadCoarseFermionEigenPack200F>("epack", ep);

    act.gauge    = "gauge";
    act.Ls       = 12;
    act.mass     = 0.000678;
    act.M5       = 1.8;
    act.b        = 1.5;
    act.c        = 0.5;
    act.boundary = "1 1 1 1";
    act.twist    = "0 0 0 0";
    application.createModule<MAction::MobiusDWF>("action", act);
}

static void addLowModes(Application &application, const std::string &name,
                        const std::string &output, const std::string &ioSplit)
{
    MUtilities::A2ALowModesCoarsePar par;

    par.eigenPack = "epack";
    par.action    = "action";
    par.output    = output;
    par.ioSplit   = ioSplit;

    application.createModule<MUtilities::A2ALowModesCoarse200>(name, par);
}

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
    bool         low     = GridCmdOptionExists(argv, argv + argc, "--low");
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
    globalPar.runId                   = std::string(low ? "split_low" : "split_write")
                                        + (read ? "_read" : "");
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
        std::vector<std::string> suffixes = low ? std::vector<std::string>{"_v", "_w"}
                                                : std::vector<std::string>{""};

        for (auto &s : suffixes)
        {
            addLoad(application, "written_full" + s, root + "/full" + s, size);
            for (auto &l : layouts)
            {
                addLoad(application, "written_split_" + l + s,
                        root + "/split_" + l + s, size);
                addCompare(application, "cmp_" + l + s, "written_full" + s,
                           "written_split_" + l + s);
            }
        }
    }
    else if (low)
    {
        addLowModeInputs(application, size);
        addLowModes(application, "low_full", root + "/full", "");
        for (auto &l : layouts)
        {
            addLowModes(application, "low_split_" + l, root + "/split_" + l,
                        layoutToPar(l));
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
