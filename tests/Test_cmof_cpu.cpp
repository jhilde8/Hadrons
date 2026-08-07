// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_cmof_cpu.cpp
 *
 * Hadrons application: runs A2AChromoMagneticOperatorField (A2ASpatialSum/
 * GEMM path) and A2AChromoMagneticOperatorFieldMT (CPU MT reference path)
 * against the same random A2A vectors and gauge field in a single run, for
 * a direct same-job correctness (and, incidentally, timing) comparison --
 * same pattern as Test_mf_cpu.cpp / Test_emf_cpu.cpp.
 *
 * Output files: cmof_gpu_out.<traj>/, cmof_mt_out.<traj>/
 * Compare with h5diff, e.g.:
 *   h5diff cmof_gpu_out.0/parity0_GitSit.h5 cmof_mt_out.0/parity0_GitSit.h5 \
 *          /parity0_GitSit /parity0_GitSit
 *
 * Usage:
 *   mpirun -n 1 ./Test_cmof_cpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
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
    // Global parameters.
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 0;
    globalPar.trajCounter.end      = 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "cmof_regression";
    globalPar.scheduler.schedulerType = "naive";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments.
    // ------------------------------------------------------------------
    int         N_i        = 8;
    int         N_j        = 8;
    int         block      = 8;
    int         cacheBlock = 8;
    std::string parities   = "0 1";
    std::string ifOrthogs  = "0 1";
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--block"))
        block      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--block"));
    if (GridCmdOptionExists(argv, argv + argc, "--cacheBlock"))
        cacheBlock = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--cacheBlock"));
    if (GridCmdOptionExists(argv, argv + argc, "--parities"))
        parities   = GridCmdOptionPayload(argv, argv + argc, "--parities");
    if (GridCmdOptionExists(argv, argv + argc, "--ifOrthogs"))
        ifOrthogs  = GridCmdOptionPayload(argv, argv + argc, "--ifOrthogs");

    // ------------------------------------------------------------------
    // Random gauge configuration and A2A vectors -- generated once,
    // referenced by name from both modules below, so both contract
    // literally identical input data in this one run.
    // ------------------------------------------------------------------
    application.createModule<MGauge::Random>("gauge");

    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // A2AChromoMagneticOperatorField -- A2ASpatialSum/GEMM path
    // ------------------------------------------------------------------
    MContraction::A2AChromoMagneticOperatorFieldPar cmofPar;
    cmofPar.block      = block;
    cmofPar.cacheBlock = cacheBlock;
    cmofPar.parities   = parities;
    cmofPar.left       = "left";
    cmofPar.right      = "right";
    cmofPar.gauge      = "gauge";
    cmofPar.output     = "cmof_gpu_out";
    cmofPar.ifOrthogs  = ifOrthogs;

    application.createModule<MContraction::A2AChromoMagneticOperatorField>("cmof_gpu", cmofPar);

    // ------------------------------------------------------------------
    // A2AChromoMagneticOperatorFieldMT -- deprecated CPU MT reference path
    // ------------------------------------------------------------------
    MContraction::A2AChromoMagneticOperatorFieldMTPar cmofMtPar;
    cmofMtPar.cacheBlock = cacheBlock;
    cmofMtPar.parities   = parities;
    cmofMtPar.left       = "left";
    cmofMtPar.right      = "right";
    cmofMtPar.gauge      = "gauge";
    cmofMtPar.output     = "cmof_mt_out";
    cmofMtPar.ifOrthogs  = ifOrthogs;

    application.createModule<MContraction::A2AChromoMagneticOperatorFieldMT>("cmof_mt", cmofMtPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
