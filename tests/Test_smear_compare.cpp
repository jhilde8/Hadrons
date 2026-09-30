// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_smear_compare.cpp
 *
 * Hadrons application: smears the same random A2A vectors with both covariant
 * smearing modules, on each of three gauge fields, and builds one meson field
 * per combination with the same contraction module. Six meson fields plus an
 * unsmeared baseline, from one set of vectors in one run.
 *
 *   modules  A2ACovariantSmearMT  Cshift path via Grid CovariantSmearing
 *            A2ACovariantSmear    6 point Laplacian stencil
 *
 *   gauge    raw    MGauge::Random, unsmeared
 *            stout  MGauge::StoutSmearing
 *            ape    MGauge::APESmear
 *
 * Comparing the two modules at fixed gauge tests the vector smearing.
 * Comparing the three gauge fields at fixed module tests what the gauge
 * smearing does to the covariant Laplacian, which is the reason the smeared
 * gauge arms are here. stout separates a specific problem with APESmear from
 * any consequence of smoothing the links at all.
 *
 * The raw field is random rather than unit because stout and APE smearing of
 * unit links is a no-op -- the staples of unit links are unit -- so all three
 * gauge arms would coincide and the test would say nothing.
 *
 * A2ACovariantSmearOrig is not included: it agrees with MT bitwise, both being
 * calls into the same Grid routine, so it adds no coverage, and its result
 * lands back in its own input array where the dependency graph cannot see it.
 *
 * The unsmeared baseline is the denominator for the elementwise ratio
 * MF_smeared/MF_unsmeared. That ratio is the quantity that exposes a smearing
 * behaving as an overall constant: a real smearing changes it element by
 * element, one that has collapsed to a scalar leaves it flat.
 *
 * Start at --N 1. One iteration is analytic,
 * MF(1) = MF(0) + coeff*<w|Gamma Laplacian|v> with coeff = alpha^2/(4N), so the
 * two modules must agree to reduction order roundoff and a disagreement
 * localises immediately. Then walk N up: a discrepancy appearing only at large
 * N points at boundary or halo handling accumulating over iterations.
 *
 * orthog_axis is fixed at 3 throughout. A2ACovariantSmear ignores the parameter
 * (its stencil hardcodes x, y, z and the -6 diagonal) while MT honours it, so
 * the two agree only at 3 and any other value would compare different
 * operators.
 *
 * Output: <output>_unsm/, and <output>_{mt,st}_{raw,stout,ape}/
 * Compare with h5diff, e.g.
 *   h5diff smear_out_mt_ape/Gamma5_0_0_0.h5 smear_out_st_ape/Gamma5_0_0_0.h5 \
 *          /Gamma5_0_0_0 /Gamma5_0_0_0
 * or take the ratio of any arm against smear_out_unsm/ to see whether the
 * smearing is acting elementwise or as a constant.
 *
 * Usage:
 *   mpirun -n 1 ./Test_smear_compare --grid 4.4.4.8 --mpi 1.1.1.1 --N 1
 */

#define HADRONS_A2AM_IO_TYPE ComplexD
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

using namespace Grid;
using namespace Hadrons;

// APE mixing weight, and the stout equivalent. Not command line options: they
// are not interchangeable with each other and tuning them is a source edit.
// The APE alpha is dimensionless and does not scale with the lattice spacing.
static const double APE_ALPHA = 0.615384615;
static const double STOUT_RHO = 0.1;

// One gauge field the vectors can be smeared with.
struct GaugeArm
{
    std::string object;
    std::string tag;
};

// One smearing module, with the parameters every arm shares. orthog_axis is
// pinned here rather than passed in; see the header note.
template <typename SmearModule, typename SmearPar>
static void addSmear(Application &application, const std::string &name,
                     const std::string &vecs, const std::string &gauge,
                     double alpha, unsigned int N)
{
    SmearPar par;

    par.a2aVectors  = vecs;
    par.gauge       = gauge;
    par.alpha       = alpha;
    par.N           = N;
    par.orthog_axis = 3;
    par.output      = "";
    par.multiFile   = false;

    application.createModule<SmearModule>(name, par);
}

// An exact copy of an A2A vector array, used to hand A2ACovariantSmear its own
// input. A2ACovariantSmearMT at alpha = 0 leaves the field untouched:
// GaussianSmear forms coeff = alpha^2/(4N) = 0 and its single iteration is
// chi = chi + 0*psi.
static void addCopy(Application &application, const std::string &name,
                    const std::string &vecs, const std::string &gauge)
{
    addSmear<MUtilities::A2ACovariantSmearMT,
             MUtilities::A2ACovariantSmearMTPar>(application, name, vecs,
                                                 gauge, 0., 1);
}

static void addMesonField(Application &application, const std::string &name,
                          const std::string &left, const std::string &right,
                          const std::string &output, int leftBlock,
                          int rightBlock, const std::string &gammas,
                          const std::vector<std::string> &mom)
{
    MContraction::A2AMesonFieldPar par;

    par.leftBlock   = leftBlock;
    par.rightBlock  = rightBlock;
    par.left        = left;
    par.right       = right;
    par.output      = output;
    par.gammas      = gammas;
    par.mom         = mom;
    par.timeSliceIO = false;

    application.createModule<MContraction::A2AMesonField>(name, par);
}

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
    // Global parameters. Naive scheduling for a deterministic module order
    // across runs; the dependency graph is fully expressed, so it is a
    // convenience rather than a correctness requirement.
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start       = 0;
    globalPar.trajCounter.end         = 1;
    globalPar.trajCounter.step        = 1;
    globalPar.runId                   = "smear_regression";
    globalPar.scheduler.schedulerType = "naive";
    globalPar.genetic.maxGen          = 1000;
    globalPar.genetic.maxCstGen       = 200;
    globalPar.genetic.popSize         = 20;
    globalPar.genetic.mutationRate    = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments.
    // ------------------------------------------------------------------
    int          N_i         = 8;
    int          N_j         = 8;
    int          leftBlock   = 8;
    int          rightBlock  = 8;
    int          momShell    = 0;
    double       alpha       = 4.1;
    unsigned int nSmear      = 1;
    // Gauge smearing steps, shared by the stout and APE arms so the two are
    // compared at equal iteration count. 45 is 48I's 25 scaled as a^-2.
    unsigned int nGauge      = 45;
    std::string  gammas      = "Gamma5";
    std::string  output_path = "smear_out";

    if (GridCmdOptionExists(argv, argv + argc, "--Ni"))
        N_i         = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Ni"));
    if (GridCmdOptionExists(argv, argv + argc, "--Nj"))
        N_j         = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nj"));
    if (GridCmdOptionExists(argv, argv + argc, "--leftBlock"))
        leftBlock   = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--leftBlock"));
    if (GridCmdOptionExists(argv, argv + argc, "--rightBlock"))
        rightBlock  = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--rightBlock"));
    if (GridCmdOptionExists(argv, argv + argc, "--mom"))
        momShell    = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--mom"));
    if (GridCmdOptionExists(argv, argv + argc, "--alpha"))
        alpha       = std::stod(GridCmdOptionPayload(argv, argv + argc, "--alpha"));
    if (GridCmdOptionExists(argv, argv + argc, "--N"))
        nSmear      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--N"));
    if (GridCmdOptionExists(argv, argv + argc, "--gaugeN"))
        nGauge      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--gaugeN"));
    if (GridCmdOptionExists(argv, argv + argc, "--gammas"))
        gammas      = cliListToPar(GridCmdOptionPayload(argv, argv + argc, "--gammas"));
    if (GridCmdOptionExists(argv, argv + argc, "--output"))
        output_path = GridCmdOptionPayload(argv, argv + argc, "--output");

    std::vector<std::string> momenta = momentumShells(momShell);

    // ------------------------------------------------------------------
    // The three gauge fields.
    // ------------------------------------------------------------------
    application.createModule<MGauge::Random>("gauge_raw");

    MGauge::StoutSmearingPar stoutPar;
    stoutPar.gauge     = "gauge_raw";
    stoutPar.steps     = nGauge;
    stoutPar.orthogDim = "3";   // a string in this module's Par, not an int
    stoutPar.rho       = STOUT_RHO;
    application.createModule<MGauge::StoutSmearing>("gauge_stout", stoutPar);

    MGauge::APESmearPar apePar;
    apePar.gauge       = "gauge_raw";
    apePar.alpha       = APE_ALPHA;
    apePar.N           = nGauge;
    apePar.orthog_axis = 3;
    application.createModule<MGauge::APESmear>("gauge_ape", apePar);

    // ------------------------------------------------------------------
    // Random A2A vectors, generated once and referenced by name from every
    // arm below, so all arms smear literally identical input data.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rvLeft, rvRight;
    rvLeft.size  = N_i; rvLeft.Ls  = 1; rvLeft.output  = ""; rvLeft.multiFile  = false;
    rvRight.size = N_j; rvRight.Ls = 1; rvRight.output = ""; rvRight.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("left",  rvLeft);
    application.createModule<MUtilities::RandomFermions>("right", rvRight);

    // ------------------------------------------------------------------
    // Unsmeared baseline, gauge independent.
    // ------------------------------------------------------------------
    addMesonField(application, "mf_unsm", "left", "right",
                  output_path + "_unsm", leftBlock, rightBlock, gammas, momenta);

    // ------------------------------------------------------------------
    // Both modules on each gauge field. MT is non-destructive and reads the
    // shared arrays; A2ACovariantSmear std::moves its input away, so it gets
    // private copies -- sharing an array with it would leave whichever arm ran
    // second smearing an empty one.
    // ------------------------------------------------------------------
    const std::vector<GaugeArm> gauges = {{"gauge_raw",   "raw"},
                                          {"gauge_stout", "stout"},
                                          {"gauge_ape",   "ape"}};

    for (auto &g : gauges)
    {
        addSmear<MUtilities::A2ACovariantSmearMT,
                 MUtilities::A2ACovariantSmearMTPar>(
            application, "l_mt_" + g.tag, "left", g.object, alpha, nSmear);
        addSmear<MUtilities::A2ACovariantSmearMT,
                 MUtilities::A2ACovariantSmearMTPar>(
            application, "r_mt_" + g.tag, "right", g.object, alpha, nSmear);
        addMesonField(application, "mf_mt_" + g.tag,
                      "l_mt_" + g.tag, "r_mt_" + g.tag,
                      output_path + "_mt_" + g.tag,
                      leftBlock, rightBlock, gammas, momenta);

        addCopy(application, "l_cp_" + g.tag, "left",  "gauge_raw");
        addCopy(application, "r_cp_" + g.tag, "right", "gauge_raw");
        addSmear<MUtilities::A2ACovariantSmear,
                 MUtilities::A2ACovariantSmearPar>(
            application, "l_st_" + g.tag, "l_cp_" + g.tag, g.object, alpha, nSmear);
        addSmear<MUtilities::A2ACovariantSmear,
                 MUtilities::A2ACovariantSmearPar>(
            application, "r_st_" + g.tag, "r_cp_" + g.tag, g.object, alpha, nSmear);
        addMesonField(application, "mf_st_" + g.tag,
                      "l_st_" + g.tag, "r_st_" + g.tag,
                      output_path + "_st_" + g.tag,
                      leftBlock, rightBlock, gammas, momenta);
    }

    application.run();

    Grid_finalize();

    return EXIT_SUCCESS;
}
