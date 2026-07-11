// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_mf_gpu.cpp
 *
 * Hadrons application: runs A2AMesonField (GPU path via A2ASpatialSum)
 * on random A2A vectors and writes HDF5 output to mf_gpu_out.<traj>/.
 *
 * Run alongside Test_fmf_cpu.cpp (same seed -> identical random vectors) then
 * diff the HDF5 outputs to validate agreement between CPU and GPU paths.
 *
 * Usage:
 *   mpirun -n 1 ./Test_mf_gpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
 *
 * The --seed argument must match Test_fmf_cpu.cpp to get identical random
 * vectors in both runs. The runId and trajectory counter are hard-coded to
 * match.
 *
 * Output files: mf_gpu_out.0/<gamma>_<mom>.h5
 *   e.g. mf_gpu_out.0/Gamma5_0_0_0.h5, mf_gpu_out.0/Identity_0_0_0.h5
 *
 * Compare against Test_fmf_cpu output with h5diff, e.g.:
 *   h5diff mf_gpu_out.0/Gamma5_0_0_0.h5 fmf_cpu_out.0/Gamma5_0_0_0.h5 \
 *          /Gamma5_0_0_0 /Gamma5_0_0_0
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
    // Global parameters - must match Test_mf_gpu.cpp exactly so that
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
    // Random A2A vectors
    // Names and sizes must match Test_mf_gpu.cpp exactly.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // A2AMesonField - GPU path (no A2AMatrixBlockComputation)
    // ------------------------------------------------------------------
    MContraction::A2AMesonFieldPar mfPar;
    mfPar.block  = block;
    mfPar.left   = "left";
    mfPar.right  = "right";
    mfPar.output = "mf_gpu_out";
    mfPar.gammas = gammas;
    mfPar.mom    = momenta;

    application.createModule<MContraction::A2AMesonField>("mf_gpu", mfPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
