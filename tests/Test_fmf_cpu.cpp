/*
 * Test_fmf_cpu.cpp
 *
 * Hadrons application: runs A2AFewMesonField (CPU path) on random A2A vectors
 * and writes HDF5 output to fmf_cpu_out.<traj>/.
 *
 * Run alongside Test_mf_gpu.cpp (same seed → identical random vectors)
 * then diff the HDF5 outputs to validate agreement between the two modules.
 *
 * Usage:
 *   mpirun -n 1 ./Test_fmf_cpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
 *
 * The --seed argument must match Test_mf_gpu.cpp to get identical random
 * vectors in both runs. The runId and trajectory counter are already
 * hard-coded to match.
 *
 * Note: A2AFewMesonField only writes output when HADRONS_A2AM_PARALLEL_IO
 * is defined at compile time.
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
    // Global parameters — must match Test_mf_gpu.cpp exactly so that
    // the RNG produces the same random vectors in both runs.
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 0;
    globalPar.trajCounter.end      = 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "mf_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments (must match Test_mf_gpu.cpp values
    // to get identical random vectors in both runs).
    // ------------------------------------------------------------------
    int N_i        = 16;
    int N_j        = 16;
    int block      = 16;
    int cacheBlock = 8;
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--block"))
        block      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--block"));
    if (GridCmdOptionExists(argv, argv + argc, "--cacheBlock"))
        cacheBlock = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--cacheBlock"));

    // ------------------------------------------------------------------
    // Random A2A vectors
    // Names and sizes must match Test_mf_gpu.cpp exactly.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // A2AFewMesonField — CPU path
    // ------------------------------------------------------------------
    MContraction::A2AFewMesonFieldPar fmfPar;
    fmfPar.cacheBlock = cacheBlock;
    fmfPar.block      = block;
    fmfPar.left       = "left";
    fmfPar.right      = "right";
    fmfPar.output     = "fmf_cpu_out";
    fmfPar.gammas     = "Gamma5 Identity";
    fmfPar.mom        = {"0 0 0", "0 0 1"};

    application.createModule<MContraction::A2AFewMesonField>("fmf_cpu", fmfPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
