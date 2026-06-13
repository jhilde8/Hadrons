/*
 * Test_emf_mt_cpu.cpp
 *
 * Hadrons application: runs A2AExtendedMesonFieldMT (CPU reference) on
 * random A2A vectors and writes HDF5 output to emf_mt_out.<traj>/.
 *
 * Run alongside Test_emf_gpu.cpp (same seed → identical random vectors)
 * then diff the HDF5 outputs to validate the GPU module.
 *
 * Usage:
 *   mpirun -n 1 ./Test_emf_mt_cpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
 *
 * The --seed argument must match Test_emf_gpu.cpp to get identical random
 * vectors in both runs. The runId and trajectory counter are already
 * hard-coded to match.
 */

#define HADRONS_A2AM_IO_TYPE ComplexD
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    Application application;

    // ------------------------------------------------------------------
    // Global parameters — must match Test_emf_gpu.cpp exactly so that
    // the RNG produces the same random vectors in both runs.
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 0;
    globalPar.trajCounter.end      = 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "emf_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional --Ni / --Nj from command line (default 8).
    // Must match Test_emf_gpu.cpp values for identical random vectors.
    // ------------------------------------------------------------------
    int N_i   = 8;
    int N_j   = 8;
    int Nloop = 8;
    int Ncb = 8;
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nloop"))
        Nloop = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nloop"));
    if (GridCmdOptionExists(argv, argv + argc, "--cacheBlock"))
        Ncb = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--cacheBlock"));

    // ------------------------------------------------------------------
    // Random A2A vectors
    // Names and sizes must match Test_emf_gpu.cpp exactly.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight, rvLoop;
    rvLeft.size  = N_i;   rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j;   rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;
    rvLoop.size  = Nloop; rvLoop.Ls  = 1; rvLoop.output  = ""; rvLoop.multiFile  = false;

    application.createModule<MUtilities::RandomFermions>("left",     rvLeft);
    application.createModule<MUtilities::RandomFermions>("right",    rvRight);
    application.createModule<MUtilities::RandomFermions>("loop_vw1", rvLoop);
    application.createModule<MUtilities::RandomFermions>("loop_vw2", rvLoop);

    // ------------------------------------------------------------------
    // A2AExtendedMesonFieldMT — CPU reference
    // ------------------------------------------------------------------
    MContraction::A2AExtendedMesonFieldMTPar emfPar;
    emfPar.cacheBlock = Ncb;
    emfPar.types      = "0 1 2 3";
    emfPar.left       = "left";
    emfPar.right      = "right";
    emfPar.loop_vw1   = "loop_vw1";
    emfPar.loop_vw2   = "loop_vw2";
    emfPar.output     = "emf_mt_out";
    emfPar.gammas1    = "GammaMU";
    emfPar.gammas2    = "GammaMU";

    application.createModule<MContraction::A2AExtendedMesonFieldMT>("emf_mt", emfPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
