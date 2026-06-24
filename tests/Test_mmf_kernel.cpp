// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constexpr constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_mmf_kernel.cpp
 *
 * Hadrons application: runs A2AMomMesonField (momentum-factored custom kernel
 * path) on random A2A vectors and writes HDF5 output to mmf_kernel_out.<traj>/.
 *
 * Run with the same --seed and --grid as Test_mf_gpu and Test_mmf_blas to get
 * identical random vectors across all three, then diff the HDF5 outputs:
 *
 *   h5diff mmf_kernel_out.0/Gamma5_0_0_0.h5 nmf_gpu_out.0/Gamma5_0_0_0.h5 \
 *          /Gamma5_0_0_0 /Gamma5_0_0_0
 *
 * Usage:
 *   mpirun -n <N> ./Test_mmf_kernel --grid 4.4.4.8 --mpi 1.1.1.1 \
 *                                    --seed "1 2 3 4" [--Ni 16] [--Nj 16] \
 *                                    [--block 16] [--mom 0] [--gammas "Gamma5 Identity"]
 */

#define HADRONS_A2AM_IO_TYPE ComplexD
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

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
    // Global parameters - must match Test_mf_gpu.cpp and Test_mmf_blas.cpp
    // exactly so that the RNG produces the same random vectors in all runs.
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
    // CLI arguments - must match Test_mf_gpu.cpp and Test_mmf_blas.cpp
    // values to get identical random vectors for diffing.
    // ------------------------------------------------------------------
    int         N_i      = 16;
    int         N_j      = 16;
    int         block    = 16;
    int         momShell = 0;
    std::string gammas   = "Gamma5 Identity";
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--block"))
        block    = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--block"));
    if (GridCmdOptionExists(argv, argv + argc, "--mom"))
        momShell = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--mom"));
    if (GridCmdOptionExists(argv, argv + argc, "--gammas"))
        gammas   = GridCmdOptionPayload(argv, argv + argc, "--gammas");

    std::vector<std::string> momenta = momentumShells(momShell);

    // ------------------------------------------------------------------
    // Random A2A vectors - module names and sizes must match Test_mf_gpu.cpp
    // and Test_mmf_blas.cpp exactly.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // A2AMomMesonField - custom kernel Phase 1 (SpinColorTrace) and Phase 2
    // ------------------------------------------------------------------
    MContraction::A2AMomMesonFieldPar mmfPar;
    mmfPar.block  = block;
    mmfPar.left   = "left";
    mmfPar.right  = "right";
    mmfPar.output = "mmf_kernel_out";
    mmfPar.gammas = gammas;
    mmfPar.mom    = momenta;

    application.createModule<MContraction::A2AMomMesonField>("mmf_kernel", mmfPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
