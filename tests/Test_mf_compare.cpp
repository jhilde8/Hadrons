// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_mf_compare.cpp
 *
 * Hadrons application: runs A2AMesonField (A2ASpatialSum ring path) and
 * A2AFewMesonField (CPU MT production path) against the same random A2A
 * vectors in a single run, for a direct CPU timing and correctness
 * comparison.
 *
 * A third slot for A2ANewMesonField is present but commented out, so the
 * driver becomes a three-way regression again as soon as a development
 * module is reintroduced under that name.
 *
 * Output files: mf_gpu_out.<traj>/, fmf_cpu_out.<traj>/
 * Compare with h5diff, e.g.:
 *   h5diff mf_gpu_out.0/Gamma5_0_0_0.h5 fmf_cpu_out.0/Gamma5_0_0_0.h5 \
 *          /Gamma5_0_0_0 /Gamma5_0_0_0
 *
 * Usage:
 *   mpirun -n 1 ./Test_mf_compare --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
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
    // Global parameters.
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
    // Parse optional CLI arguments.
    // ------------------------------------------------------------------
    int         N_i        = 16;
    int         N_j        = 16;
    int         block      = 16;
    int         cacheBlock = 8;
    int         momShell   = 0;
    std::string gammas     = "Gamma5 Identity";
    std::string output_path= "mf_gpu_out";
    // Presence-only flag. The oracle (A2AFewMesonField) always writes the
    // merged layout, so with this set the comparison is per-timeslice files
    // against one full-nt file -- use diffs/diff_mf_ts.py rather than
    // diff_mf.py. Needs P_t > 1 in --mpi or the flag is a no-op.
    bool        tsIO       = false;
    if (GridCmdOptionExists(argv, argv + argc, "--timeSliceIO"))
        tsIO       = true;
    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j        = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--block"))
        block      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--block"));
    if (GridCmdOptionExists(argv, argv + argc, "--cacheBlock"))
        cacheBlock = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--cacheBlock"));
    if (GridCmdOptionExists(argv, argv + argc, "--mom"))
        momShell   = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--mom"));
    if (GridCmdOptionExists(argv, argv + argc, "--gammas"))
        gammas     = GridCmdOptionPayload(argv, argv + argc, "--gammas");
    if (GridCmdOptionExists(argv, argv + argc, "--output"))
        output_path     = GridCmdOptionPayload(argv, argv + argc, "--output");

    std::vector<std::string> momenta = momentumShells(momShell);

    // ------------------------------------------------------------------
    // Random A2A vectors -- generated once, referenced by name from all
    // three modules below, so all three contract literally identical
    // input data in this one run.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // A2AMesonField -- current GPU path (A2ASpatialSum / batched GEMM)
    // ------------------------------------------------------------------
    MContraction::A2AMesonFieldPar mfPar;
    mfPar.block      = block;
    mfPar.cacheBlock = cacheBlock;
    mfPar.left   = "left";
    mfPar.right  = "right";
    mfPar.output = output_path;
    mfPar.gammas = gammas;
    mfPar.mom    = momenta;
    mfPar.timeSliceIO = tsIO;

    application.createModule<MContraction::A2AMesonField>("mf_gpu", mfPar);

    // ------------------------------------------------------------------
    // A2AFewMesonField -- CPU MT production path (contract-once + mac)
    // ------------------------------------------------------------------
    MContraction::A2AFewMesonFieldPar fmfPar;
    fmfPar.cacheBlock = cacheBlock;
    fmfPar.block      = block;
    fmfPar.left       = "left";
    fmfPar.right      = "right";
    fmfPar.output     = "fmf_cpu_out";
    fmfPar.gammas     = gammas;
    fmfPar.mom        = momenta;

    application.createModule<MContraction::A2AFewMesonField>("fmf_cpu", fmfPar);

    // ------------------------------------------------------------------
    // A2ANewMesonField -- the development slot. No such module exists right
    // now; reintroduce one under that name to compare a candidate against
    // both A2AMesonField and A2AFewMesonField in a single run, and
    // uncomment this block. Left commented out so the driver runs as the
    // two-way FMF vs MF comparison in the meantime.
    // ------------------------------------------------------------------
    // MContraction::A2ANewMesonFieldPar nmfPar;
    // nmfPar.block      = block;
    // nmfPar.cacheBlock = cacheBlock;
    // nmfPar.left   = "left";
    // nmfPar.right  = "right";
    // nmfPar.output = "mf_new_out";
    // nmfPar.gammas = gammas;
    // nmfPar.mom    = momenta;
    //
    // application.createModule<MContraction::A2ANewMesonField>("mf_new", nmfPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
