// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_smear_compare.cpp
 *
 * Hadrons application: runs the A2A covariant smearing modules over the same
 * random A2A vectors and the same gauge field in a single run, then builds one
 * meson field per smearing using the SAME contraction module. The smearing is
 * the only variable, so any difference between the output meson fields is a
 * difference between the smearing implementations.
 *
 *   A2ACovariantSmearMT    - Cshift path via Grid CovariantSmearing, production
 *   A2ACovariantSmear      - 6 point Laplacian stencil, newer
 *   A2ACovariantSmearOrig  - Cshift path, smears in place, production (--orig)
 *
 * An unsmeared meson field is written too. It is the reference for the
 * elementwise ratio MF_smeared/MF_unsmeared, and that ratio is the quantity
 * that exposes a smearing behaving as an overall constant: a real smearing
 * changes it element by element, one that has collapsed to a scalar leaves it
 * flat.
 *
 * GAUGE CHOICE - run both. With --gauge unit the covariant Laplacian reduces
 * to the free one, whose eigenvectors are plane waves, so the zero momentum
 * component is an exact eigenvector with smearing eigenvalue 1 and a correct
 * module must preserve it exactly. That makes unit links the sharper test of
 * the shift and halo logic, but every link is the identity, so a wrong gauge
 * index or a mis-multiplied link is invisible. --gauge random exercises the
 * parallel transport that unit links hide.
 *
 * ITERATION COUNT - start at --N 1. One iteration is analytic,
 * MF(1) = MF(0) + coeff*<w|Gamma Laplacian|v> with coeff = alpha^2/(4N), so
 * all modules must agree to reduction order roundoff and a disagreement
 * localises immediately. Then walk N up: a discrepancy that only appears at
 * large N points at boundary or halo handling accumulating over iterations
 * rather than at the arithmetic of a single application.
 *
 * orthog_axis is fixed at 3 and is deliberately not a CLI knob.
 * A2ACovariantSmear ignores the parameter (its stencil hardcodes x, y, z and
 * the -6 diagonal) while MT and Orig honour it, so the three agree only at 3
 * and any other value would compare different operators.
 *
 * Output: <output>_unsm/, <output>_mt/, <output>_stencil/, <output>_orig/
 * Compare with h5diff, e.g.
 *   h5diff smear_out_mt/Gamma5_0_0_0.h5 smear_out_stencil/Gamma5_0_0_0.h5 \
 *          /Gamma5_0_0_0 /Gamma5_0_0_0
 * or take the ratio against smear_out_unsm/ to see whether the smearing is
 * acting elementwise or as a constant.
 *
 * Usage:
 *   mpirun -n 1 ./Test_smear_compare --grid 4.4.4.8 --mpi 1.1.1.1 \
 *          --alpha 4.1 --N 1 --gauge unit
 */

#define HADRONS_A2AM_IO_TYPE ComplexD
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

using namespace Grid;
using namespace Hadrons;

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

// An exact copy of an A2A vector array, used to hand the destructive modules
// their own input. A2ACovariantSmearMT at alpha = 0 leaves the field untouched:
// GaussianSmear forms coeff = alpha^2/(4N) = 0 and its iterations reduce to
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
    // Global parameters. The scheduler is naive rather than genetic on
    // purpose: A2ACovariantSmearOrig writes its result back into its input
    // array and advertises only a scratch field as its output, so nothing in
    // the dependency graph forces it to run before the meson field that reads
    // that array. Insertion order is what makes the --orig arm correct.
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
    std::string  gammas      = "Gamma5";
    std::string  gaugeKind   = "unit";
    std::string  output_path = "smear_out";
    // A2ACovariantSmearOrig is opt-in: see the scheduler note above.
    bool         doOrig      = false;

    if (GridCmdOptionExists(argv, argv + argc, "--orig"))
        doOrig      = true;
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
    if (GridCmdOptionExists(argv, argv + argc, "--gammas"))
        gammas      = cliListToPar(GridCmdOptionPayload(argv, argv + argc, "--gammas"));
    if (GridCmdOptionExists(argv, argv + argc, "--gauge"))
        gaugeKind   = GridCmdOptionPayload(argv, argv + argc, "--gauge");
    if (GridCmdOptionExists(argv, argv + argc, "--output"))
        output_path = GridCmdOptionPayload(argv, argv + argc, "--output");

    std::vector<std::string> momenta = momentumShells(momShell);

    // ------------------------------------------------------------------
    // Gauge field, shared by every smearing module.
    // ------------------------------------------------------------------
    if (gaugeKind == "unit")
    {
        application.createModule<MGauge::Unit>("gauge");
    }
    else if (gaugeKind == "random")
    {
        application.createModule<MGauge::Random>("gauge");
    }
    else
    {
        HADRONS_ERROR(Argument, "--gauge must be 'unit' or 'random', got '"
                                + gaugeKind + "'");
    }

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
    // Unsmeared baseline. Also the denominator for the elementwise ratio.
    // ------------------------------------------------------------------
    addMesonField(application, "mf_unsm", "left", "right",
                  output_path + "_unsm", leftBlock, rightBlock, gammas, momenta);

    // ------------------------------------------------------------------
    // A2ACovariantSmearMT. Non-destructive, so it reads the shared arrays.
    // ------------------------------------------------------------------
    addSmear<MUtilities::A2ACovariantSmearMT,
             MUtilities::A2ACovariantSmearMTPar>(application, "l_mt", "left",
                                                 "gauge", alpha, nSmear);
    addSmear<MUtilities::A2ACovariantSmearMT,
             MUtilities::A2ACovariantSmearMTPar>(application, "r_mt", "right",
                                                 "gauge", alpha, nSmear);
    addMesonField(application, "mf_mt", "l_mt", "r_mt",
                  output_path + "_mt", leftBlock, rightBlock, gammas, momenta);

    // ------------------------------------------------------------------
    // A2ACovariantSmear. std::move empties its input, so it gets private
    // copies -- sharing "left"/"right" with it would leave whichever arm ran
    // second smearing an empty array.
    // ------------------------------------------------------------------
    addCopy(application, "l_cp", "left",  "gauge");
    addCopy(application, "r_cp", "right", "gauge");
    addSmear<MUtilities::A2ACovariantSmear,
             MUtilities::A2ACovariantSmearPar>(application, "l_st", "l_cp",
                                               "gauge", alpha, nSmear);
    addSmear<MUtilities::A2ACovariantSmear,
             MUtilities::A2ACovariantSmearPar>(application, "r_st", "r_cp",
                                               "gauge", alpha, nSmear);
    addMesonField(application, "mf_stencil", "l_st", "r_st",
                  output_path + "_stencil", leftBlock, rightBlock, gammas, momenta);

    // ------------------------------------------------------------------
    // A2ACovariantSmearOrig. Smears in place, so the meson field reads the
    // copies the module was handed rather than the module's own name.
    // ------------------------------------------------------------------
    if (doOrig)
    {
        addCopy(application, "l_cpo", "left",  "gauge");
        addCopy(application, "r_cpo", "right", "gauge");
        addSmear<MUtilities::A2ACovariantSmearOrig,
                 MUtilities::A2ACovariantSmearOrigPar>(application, "l_or",
                                                       "l_cpo", "gauge",
                                                       alpha, nSmear);
        addSmear<MUtilities::A2ACovariantSmearOrig,
                 MUtilities::A2ACovariantSmearOrigPar>(application, "r_or",
                                                       "r_cpo", "gauge",
                                                       alpha, nSmear);
        addMesonField(application, "mf_orig", "l_cpo", "r_cpo",
                      output_path + "_orig", leftBlock, rightBlock, gammas, momenta);
    }

    application.run();

    Grid_finalize();

    return EXIT_SUCCESS;
}
