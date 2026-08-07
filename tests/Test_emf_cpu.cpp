// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_emf_cpu.cpp
 *
 * Hadrons application: runs A2AExtendedMesonField (A2ASpatialSum/GEMM path)
 * and A2AExtendedMesonFieldMT (CPU MT reference path) against the same
 * random A2A vectors and quark loop in a single run, for a direct
 * same-job correctness (and, incidentally, timing) comparison -- same
 * pattern as Test_mf_cpu.cpp.
 *
 * Output files: emf_gpu_out.<traj>/, emf_mt_out.<traj>/
 * Compare with h5diff, e.g.:
 *   h5diff emf_gpu_out.0/type0_GammaMU_GammaMU.h5 emf_mt_out.0/type0_GammaMU_GammaMU.h5 \
 *          /type0_GammaMU_GammaMU /type0_GammaMU_GammaMU
 *
 * Usage:
 *   mpirun -n 1 ./Test_emf_cpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
 */

#define HADRONS_A2AM_IO_TYPE ComplexD
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Hadrons/Modules/MContraction/A2AExtendedMesonFieldMT.hpp>

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
    globalPar.runId                = "emf_regression";
    // Sequential (non-genetic) scheduler -- sidesteps GeneticScheduler
    // entirely, same reasoning as Test_emf_gpu.cpp.
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
    int         Nloop      = 4;
    int         block      = 8;
    int         cacheBlock = 8;
    std::string types      = "0 1 2 3";
    std::string gammas1    = "GammaMU";
    std::string gammas2    = "GammaMU";
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nloop"))
        Nloop      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nloop"));
    if (GridCmdOptionExists(argv, argv + argc, "--block"))
        block      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--block"));
    if (GridCmdOptionExists(argv, argv + argc, "--cacheBlock"))
        cacheBlock = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--cacheBlock"));
    if (GridCmdOptionExists(argv, argv + argc, "--types"))
        types      = GridCmdOptionPayload(argv, argv + argc, "--types");
    if (GridCmdOptionExists(argv, argv + argc, "--gammas1"))
        gammas1    = GridCmdOptionPayload(argv, argv + argc, "--gammas1");
    if (GridCmdOptionExists(argv, argv + argc, "--gammas2"))
        gammas2    = GridCmdOptionPayload(argv, argv + argc, "--gammas2");

    // ------------------------------------------------------------------
    // Random A2A vectors and quark loop -- generated once, referenced by
    // name from both modules below, so both contract literally identical
    // input data in this one run.
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
    // A2AExtendedMesonField -- A2ASpatialSum/GEMM path
    // ------------------------------------------------------------------
    MContraction::A2AExtendedMesonFieldPar emfPar;
    emfPar.block      = block;
    emfPar.cacheBlock = cacheBlock;
    emfPar.types      = types;
    emfPar.left       = "left";
    emfPar.right      = "right";
    emfPar.loop_vw1   = "loop_vw1";
    emfPar.loop_vw2   = "loop_vw2";
    emfPar.output     = "emf_gpu_out";
    emfPar.gammas1    = gammas1;
    emfPar.gammas2    = gammas2;

    application.createModule<MContraction::A2AExtendedMesonField>("emf_gpu", emfPar);

    // ------------------------------------------------------------------
    // A2AExtendedMesonFieldMT -- CPU MT reference path
    // ------------------------------------------------------------------
    MContraction::A2AExtendedMesonFieldMTPar emfMtPar;
    emfMtPar.cacheBlock = cacheBlock;
    emfMtPar.types      = types;
    emfMtPar.left       = "left";
    emfMtPar.right      = "right";
    emfMtPar.loop_vw1   = "loop_vw1";
    emfMtPar.loop_vw2   = "loop_vw2";
    emfMtPar.output     = "emf_mt_out";
    emfMtPar.gammas1    = gammas1;
    emfMtPar.gammas2    = gammas2;

    application.createModule<MContraction::A2AExtendedMesonFieldMT>("emf_mt", emfMtPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
