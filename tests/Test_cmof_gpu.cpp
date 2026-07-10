// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_cmof_gpu.cpp
 *
 * Hadrons application: runs A2AChromoMagneticOperatorField (GPU path) on
 * random A2A vectors with a random gauge field and writes HDF5 output to
 * cmof_gpu_out.<traj>/.
 *
 * Run alongside Test_cmof_mt_cpu.cpp (same seed -> identical random vectors
 * and gauge configuration) then diff the HDF5 outputs to validate the GPU
 * module.
 *
 * Usage:
 *   mpirun -n 1 ./Test_cmof_gpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
 *
 * The --seed argument must match Test_cmof_mt_cpu.cpp to get identical random
 * fields in both runs. The runId and trajectory counter are already
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
    // Global parameters - must match Test_cmof_mt_cpu.cpp exactly so that
    // the RNG produces the same random fields in both runs.
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 0;
    globalPar.trajCounter.end      = 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "cmof_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments (must match Test_cmof_mt_cpu.cpp values
    // to get identical random fields in both runs).
    // ------------------------------------------------------------------
    int N_i        = 8;
    int N_j        = 8;
    int cacheBlock = 3;
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--cacheBlock"))
        cacheBlock = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--cacheBlock"));

    // ------------------------------------------------------------------
    // Random gauge configuration
    // Name must match Test_cmof_mt_cpu.cpp exactly.
    // ------------------------------------------------------------------
    application.createModule<MGauge::Random>("gauge");

    // ------------------------------------------------------------------
    // Random A2A vectors
    // Names and sizes must match Test_cmof_mt_cpu.cpp exactly.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // A2AChromoMagneticOperatorField - GPU path
    // ------------------------------------------------------------------
    MContraction::A2AChromoMagneticOperatorFieldPar cmofPar;
    cmofPar.cacheBlock = cacheBlock;
    cmofPar.parities   = "0 1";
    cmofPar.left       = "left";
    cmofPar.right      = "right";
    cmofPar.gauge      = "gauge";
    cmofPar.output     = "cmof_gpu_out";
    cmofPar.ifOrthogs  = "0 1";

    application.createModule<MContraction::A2AChromoMagneticOperatorField>("cmof_gpu", cmofPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
