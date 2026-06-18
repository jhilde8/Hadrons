// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_covsmear.cpp
 *
 * Hadrons application: runs A2ACovariantSmear (stencil path) and
 * A2ACovariantSmearMT (sequential Cshift/GaussianSmear path) on the same
 * random A2A vectors and random gauge field, then contracts each smeared set
 * into a meson field via A2ANewMesonField and writes HDF5 output.
 *
 * h5diff the two meson field outputs to verify that both smearing paths
 * produce identical results.
 *
 * Usage:
 *   mpirun -n 1 ./Test_covsmear --grid 4.4.4.8 --mpi 1.1.1.1
 *
 * Output directories:
 *   covsmear_stencil_mf_out.<traj>/   (stencil smearing -> meson field)
 *   covsmear_mt_mf_out.<traj>/        (Cshift smearing  -> meson field)
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

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 0;
    globalPar.trajCounter.end      = 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "covsmear_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    int         Ni     = 4;
    int         N      = 3;
    double      alpha  = 2.0;
    int         block  = 4;
    std::string gammas = "Gamma5";

    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        Ni    = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--N"))
        N     = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--N"));
    if (GridCmdOptionExists(argv, argv + argc, "--alpha"))
        alpha = std::stod(GridCmdOptionPayload(argv, argv + argc, "--alpha"));
    if (GridCmdOptionExists(argv, argv + argc, "--block"))
        block = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--block"));
    if (GridCmdOptionExists(argv, argv + argc, "--gammas"))
        gammas = GridCmdOptionPayload(argv, argv + argc, "--gammas");

    // ------------------------------------------------------------------
    // Random gauge configuration
    // ------------------------------------------------------------------
    application.createModule<MGauge::Random>("gauge");

    // ------------------------------------------------------------------
    // Random A2A vectors (shared input for both smearing paths)
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvPar;
    rvPar.size      = Ni;
    rvPar.Ls        = 1;
    rvPar.output    = "";
    rvPar.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("vecs", rvPar);

    // ------------------------------------------------------------------
    // A2ACovariantSmear - stencil path (module under test)
    // ------------------------------------------------------------------
    MUtilities::A2ACovariantSmearPar stencilPar;
    stencilPar.a2aVectors  = "vecs";
    stencilPar.gauge       = "gauge";
    stencilPar.alpha       = alpha;
    stencilPar.N           = N;
    stencilPar.orthog_axis = Nd - 1;
    stencilPar.output      = "";
    stencilPar.multiFile   = false;

    application.createModule<MUtilities::A2ACovariantSmear>("smear_stencil", stencilPar);

    // ------------------------------------------------------------------
    // A2ACovariantSmearMT - sequential Cshift reference path
    // ------------------------------------------------------------------
    MUtilities::A2ACovariantSmearMTPar mtPar;
    mtPar.a2aVectors  = "vecs";
    mtPar.gauge       = "gauge";
    mtPar.alpha       = alpha;
    mtPar.N           = N;
    mtPar.orthog_axis = Nd - 1;
    mtPar.output      = "";
    mtPar.multiFile   = false;

    application.createModule<MUtilities::A2ACovariantSmearMT>("smear_mt", mtPar);

    // ------------------------------------------------------------------
    // A2ANewMesonField on stencil-smeared vectors
    // ------------------------------------------------------------------
    MContraction::A2ANewMesonFieldPar mfStencilPar;
    mfStencilPar.block  = block;
    mfStencilPar.left   = "smear_stencil";
    mfStencilPar.right  = "smear_stencil";
    mfStencilPar.output = "covsmear_stencil_mf_out";
    mfStencilPar.gammas = gammas;
    mfStencilPar.mom    = {"0 0 0"};

    application.createModule<MContraction::A2ANewMesonField>("mf_stencil", mfStencilPar);

    // ------------------------------------------------------------------
    // A2ANewMesonField on MT-smeared vectors
    // ------------------------------------------------------------------
    MContraction::A2ANewMesonFieldPar mfMtPar;
    mfMtPar.block  = block;
    mfMtPar.left   = "smear_mt";
    mfMtPar.right  = "smear_mt";
    mfMtPar.output = "covsmear_mt_mf_out";
    mfMtPar.gammas = gammas;
    mfMtPar.mom    = {"0 0 0"};

    application.createModule<MContraction::A2ANewMesonField>("mf_mt", mfMtPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
