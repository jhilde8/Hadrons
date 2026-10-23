// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_split_load.cpp
 *
 * Hadrons application: MIO::LoadCombinedA2AVecsV with ioSplit against the
 * same module without it (the full-grid read) on the 64I light V vectors
 * already on disk, one arm per --io-split layout. The file stems are the production ones from
 * PyXmlGen/config.py (low_filestem("l", "v"), high_stem("l")); edit them here
 * if the tree moves.
 *
 * Correctness (default). Every arm reads the same files -- --n-low low modes,
 * then --high-size from each --high-ext -- and MUtilities::CompareA2AVecs
 * checks each split arm against the full-grid one. Pass is exact equality:
 * both paths read the same bytes, Grid_unsplit only copies them, and the
 * 1/nHit multiply is the same on both sides. The defaults (20 low, 24 per
 * extension, two extensions) are not multiples of any partition count up to
 * 16, so the last batch of each block has idle partitions, and the second
 * extension exercises the per-block offsets.
 *
 * Timing (--timing). The full-grid arm reads --high-ext entry 0, then one
 * check arm per layout reads the same entry and is compared against it, so
 * the expensive full-grid read is shared with the correctness check. The
 * check arms read warm files, so their timings mean nothing. The timed split
 * arms come last, layout k on entry k + 1 alone, so none reads a file another
 * arm has just pulled into the server caches. No low modes. Needs one
 * extension per layout plus one; the default is all eight sloppy light hits,
 * so up to seven layouts. Read the per-module timers: "read", "read wait"
 * and "unsplit" on the v_split_* arms, "read" on v_full.
 *
 * Layouts are dot separated, several comma separated. Each must divide --mpi
 * dimension by dimension.
 *
 * Usage:
 *   srun ... ./Test_split_load --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --io-split 2.2.2.4
 *   srun ... ./Test_split_load --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --timing --high-size 256 --io-split 2.2.2.4,2.2.4.4,4.4.4.4
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

using namespace Grid;
using namespace Hadrons;

static const std::string LOW_STEM  = "/lustre/orion/phy157/proj-shared/phy157_dwf/jhilde/main_64I/vw/low/lo_v";
static const std::string HIGH_STEM = "/lustre/orion/phy157/proj-shared/phy157_dwf/jhilde/main_64I/vw/sloppy/";

// What every arm shares; each arm picks its own extensions and nLow.
struct LoadSpec
{
    unsigned int highSize;
    unsigned int nHit;
};

static void addFullLoad(Application &application, const std::string &name,
                        const LoadSpec &spec, const unsigned int nLow,
                        const std::vector<std::string> &ext)
{
    MIO::LoadCombinedA2AVecsVPar par;

    par.lowFilestem    = LOW_STEM;
    par.nLow           = nLow;
    par.highStem       = HIGH_STEM;
    par.highExtensions = ext;
    par.highSize       = spec.highSize;
    par.nHit           = spec.nHit;

    application.createModule<MIO::LoadCombinedA2AVecsV>(name, par);
}

static void addSplitLoad(Application &application, const std::string &name,
                         const LoadSpec &spec, const unsigned int nLow,
                         const std::vector<std::string> &ext,
                         const std::string &ioSplit)
{
    MIO::LoadCombinedA2AVecsVPar par;

    par.lowFilestem    = LOW_STEM;
    par.nLow           = nLow;
    par.highStem       = HIGH_STEM;
    par.highExtensions = ext;
    par.highSize       = spec.highSize;
    par.nHit           = spec.nHit;
    par.ioSplit        = ioSplit;

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

// "2.2.2.4" -> "2 2 2 4", the form strToVec reads.
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
    int          traj      = 1540;
    unsigned int nLow      = 20;
    bool         timing    = GridCmdOptionExists(argv, argv + argc, "--timing");
    std::string  extList   = timing ? "l0_v,l1_v,l2_v,l3_v,l4_v,l5_v,l6_v,l7_v"
                                    : "l0_v,l1_v";
    std::string  splitList;
    LoadSpec     spec{24, 8};

    if (GridCmdOptionExists(argv, argv + argc, "--traj"))
        traj          = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--traj"));
    if (GridCmdOptionExists(argv, argv + argc, "--n-low"))
        nLow          = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--n-low"));
    if (GridCmdOptionExists(argv, argv + argc, "--high-size"))
        spec.highSize = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--high-size"));
    if (GridCmdOptionExists(argv, argv + argc, "--n-hit"))
        spec.nHit     = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--n-hit"));
    if (GridCmdOptionExists(argv, argv + argc, "--high-ext"))
        extList       = GridCmdOptionPayload(argv, argv + argc, "--high-ext");
    if (GridCmdOptionExists(argv, argv + argc, "--io-split"))
        splitList     = GridCmdOptionPayload(argv, argv + argc, "--io-split");

    std::vector<std::string> ext     = strToVec<std::string>(cliListToPar(extList));
    std::vector<std::string> layouts = strToVec<std::string>(cliListToPar(splitList));

    if (layouts.empty())
    {
        HADRONS_ERROR(Argument, "--io-split needs at least one layout");
    }
    if (timing && ext.size() < layouts.size() + 1)
    {
        HADRONS_ERROR(Argument, "--timing needs one --high-ext per arm ("
                      + std::to_string(layouts.size() + 1) + ")");
    }

    // ------------------------------------------------------------------
    // Global parameters. Naive scheduling keeps the arms in creation order,
    // which the timing mode relies on to put the full-grid arm first.
    // ------------------------------------------------------------------
    Application application;
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start       = traj;
    globalPar.trajCounter.end         = traj + 1;
    globalPar.trajCounter.step        = 1;
    globalPar.runId                   = "split_load";
    globalPar.scheduler.schedulerType = "naive";
    globalPar.genetic.maxGen          = 1000;
    globalPar.genetic.maxCstGen       = 200;
    globalPar.genetic.popSize         = 20;
    globalPar.genetic.mutationRate    = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // The arms.
    // ------------------------------------------------------------------
    if (timing)
    {
        addFullLoad(application, "v_full", spec, 0, {ext[0]});
        for (auto &l : layouts)
        {
            addSplitLoad(application, "v_check_" + l, spec, 0, {ext[0]},
                         layoutToPar(l));
            addCompare(application, "cmp_" + l, "v_full", "v_check_" + l);
        }
        for (unsigned int k = 0; k < layouts.size(); ++k)
        {
            addSplitLoad(application, "v_split_" + layouts[k], spec, 0,
                         {ext[k + 1]}, layoutToPar(layouts[k]));
        }
    }
    else
    {
        addFullLoad(application, "v_full", spec, nLow, ext);
        for (auto &l : layouts)
        {
            addSplitLoad(application, "v_split_" + l, spec, nLow, ext,
                         layoutToPar(l));
            addCompare(application, "cmp_" + l, "v_full", "v_split_" + l);
        }
    }

    application.run();

    Grid_finalize();

    return EXIT_SUCCESS;
}
