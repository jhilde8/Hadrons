// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_a2a_evalfpe_gpu.cpp
 *
 * Reproduces the MSolver::A2AVectors module pattern used for the light quark
 * in production (noise_ud -> a2a_l, mixed-precision solver mcg_l,
 * eigenPack epackD) but with a random gauge field and a synthetic, in-memory
 * eigenpack instead of loading the real ~2200-file eigenpack from disk. This
 * lets makeLowModeV/makeLowModeW's eigenvalue division be tested directly
 * and cheaply, without paying the multi-minute production IO cost.
 *
 * Action/mass here use the standard 64I ensemble parameters (Ls=12, M5=1.8,
 * b=1.5, c=0.5, mass=0.000678) with standard (non-z) Mobius DWF -- adjust
 * b/c below if your local 64I config differs; they only affect solver
 * conditioning, not the eval-division math under test.
 *
 * mcg_l mirrors production's mixed-precision solver: a single-precision
 * inner action (mdwf_lf, on a single-precision cast of the gauge field) and
 * a double-precision outer action (mdwf_l), per
 * Hadrons/Modules/MSolver/MixedPrecisionRBPrecCG.hpp. Its outerGuesser is
 * wired to an MGuesser::ExactDeflation guesser built from the same epackD
 * eigenpack, matching production (which always deflates the high-mode
 * solve). Leaving outerGuesser empty falls back to Grid's ZeroGuesser,
 * which -- unlike DeflatedGuesser/ExactDeflatedGuesser -- never sets
 * guess.Checkerboard(); combined with A2AVectors' automatic use of the
 * "_subtract" solver whenever Nl_ > 0, that mistagged zero field survives
 * into SchurRedBlackBase's subGuess subtraction and trips Lattice_ET.h's
 * checkerboard-consistency GRID_ASSERT. That's a real latent Grid bug but
 * unrelated to the eval-division FPE this test targets, so it's avoided
 * here by deflating for real rather than patching Grid's ZeroGuesser.
 * makeLowModeV/W (what's under test) still runs first in
 * TA2AVectors::execute(), before the high-mode solve. If the eval division
 * is the culprit this should crash before the solver is ever invoked. If
 * instead you hit a "CG did not converge" HADRONS_ERROR, that's the
 * mixed-precision CG struggling on a non-thermalized random gauge field --
 * a normal exception, not a signal, and unrelated to the FPE hypothesis;
 * loosen --residual or raise maxInnerIteration/maxOuterIteration below if
 * it comes up before you've ruled out the low-mode path.
 *
 * The eigenvalues are spread linearly between --evalMin and --evalMax, with
 * index 0 holding the smallest (most dangerous) value. Lower --evalMin
 * (down to and including 0.0) to probe whether/where the 1/eval division in
 * A2AVectors.hpp trips a floating point exception.
 *
 * Run with --debug-signals so Grid installs the SIGFPE trap/backtrace
 * handler (same flag used in production).
 *
 * Usage (full 64I volume, needs a real multi-node allocation):
 *   srun ... ./Test_a2a_evalfpe_gpu --grid 64.64.64.128 --mpi <MxNxPxQ> \
 *       --seed "1 2 3 4" --debug-signals --Nl 4 --evalMin 1e-4 --evalMax 2.0
 *
 * For a quick single-rank smoke test before running at full scale:
 *   mpirun -n 1 ./Test_a2a_evalfpe_gpu --grid 4.4.4.4 --mpi 1.1.1.1 \
 *       --seed "1 2 3 4" --debug-signals --Nl 4 --evalMin 1e-4 --evalMax 2.0
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// ---------------------------------------------------------------------------
// Minimal in-memory replacement for MIO::LoadFermionEigenPack: fills a
// BaseFermionEigenPack<FImpl> with random eigenvectors and a controllable
// eigenvalue spread, instead of reading v<k>.bin files from disk.
// ---------------------------------------------------------------------------
BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MUtilities)

class RandomFermionEigenPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomFermionEigenPackPar,
                                    unsigned int, size,
                                    unsigned int, Ls,
                                    double,       evalMin,
                                    double,       evalMax);
};

template <typename FImpl>
class TRandomFermionEigenPack: public Module<RandomFermionEigenPackPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef BaseFermionEigenPack<FImpl> Pack;
public:
    TRandomFermionEigenPack(const std::string name)
    : Module<RandomFermionEigenPackPar>(name)
    {}
    virtual ~TRandomFermionEigenPack(void) {};
    virtual std::vector<std::string> getInput(void)
    {
        return {};
    }
    virtual std::vector<std::string> getOutput(void)
    {
        return {getName()};
    }
    virtual void setup(void)
    {
        envCreate(Pack, getName(), par().Ls, par().size,
                  envGetRbGrid(FermionField, par().Ls));
    }
    virtual void execute(void)
    {
        auto &epack = envGet(Pack, getName());
        // GridParallelRNG::fill's automatic checkerboard handling builds its
        // internal temporary on rng4d()'s own (4D) grid, which only matches
        // a checkerboarded *4D* target. For a 5D redblack eigenvector grid
        // that temporary is the wrong dimensionality, so fill it on the full
        // 5D grid ourselves (same pattern as RandomVectors.hpp for Ls>1) and
        // project down, exactly as makeLowModeV does internally.
        FermionField full(envGetGrid(FermionField, par().Ls));

        LOG(Message) << "Generating " << par().size
                     << " random low-mode eigenvectors, eval in ["
                     << par().evalMin << ", " << par().evalMax << "]"
                     << std::endl;
        for (unsigned int i = 0; i < par().size; ++i)
        {
            random(rng4d(), full);
            pickCheckerboard(Odd, epack.evec[i], full);
            double frac  = (par().size > 1) ?
                           double(i) / double(par().size - 1) : 0.;
            epack.eval[i] = par().evalMin + frac * (par().evalMax - par().evalMin);
            LOG(Message) << "  eval[" << i << "] = " << epack.eval[i] << std::endl;
        }
    }
};

// Plain MODULE_REGISTER (not MODULE_REGISTER_TMP): the _TMP variant emits an
// "extern template class" declaration that assumes an explicit out-of-line
// instantiation exists in a companion .cpp compiled into libHadrons (as real
// modules have). Since TRandomFermionEigenPack only exists in this test
// file, that promise goes unmet and the linker can't find the vtable. Plain
// MODULE_REGISTER lets the compiler instantiate it normally, right here.
MODULE_REGISTER(RandomFermionEigenPack, ARG(TRandomFermionEigenPack<FIMPL>), MUtilities);

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE

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
    // Global parameters
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start        = 0;
    globalPar.trajCounter.end          = 1;
    globalPar.trajCounter.step         = 1;
    globalPar.runId                    = "a2a_evalfpe_regression";
    // Modules below are created in dependency order, so the naive
    // (sequential, non-genetic) scheduler is valid here and sidesteps
    // GeneticScheduler entirely -- useful while GeneticScheduler trips
    // the --debug-signals FPE trap (see doMutation()/selectPair()).
    globalPar.scheduler.schedulerType  = "naive";
    globalPar.genetic.maxGen           = 1000;
    globalPar.genetic.maxCstGen        = 200;
    globalPar.genetic.popSize          = 20;
    globalPar.genetic.mutationRate     = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments
    // ------------------------------------------------------------------
    unsigned int Nl      = 4;
    double       evalMin = 1e-4;
    double       evalMax = 2.0;
    if (GridCmdOptionExists(argv, argv + argc, "--Nl"))
        Nl = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nl"));
    if (GridCmdOptionExists(argv, argv + argc, "--evalMin"))
        evalMin = std::stod(GridCmdOptionPayload(argv, argv + argc, "--evalMin"));
    if (GridCmdOptionExists(argv, argv + argc, "--evalMax"))
        evalMax = std::stod(GridCmdOptionPayload(argv, argv + argc, "--evalMax"));

    const unsigned int Ls = 12;

    // ------------------------------------------------------------------
    // Random gauge field (stand-in for MIO::LoadNersc), plus single
    // precision cast for the mixed-precision solver's inner action.
    // ------------------------------------------------------------------
    application.createModule<MGauge::Random>("gauge");

    MUtilities::PrecisionCastPar castPar;
    castPar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugef", castPar);

    // ------------------------------------------------------------------
    // Mobius DWF action, standard 64I parameters -- double precision
    // (outer) and single precision (inner) copies for the mixed-precision
    // solver below.
    // ------------------------------------------------------------------
    MAction::MobiusDWFPar actionPar;
    actionPar.gauge    = "gauge";
    actionPar.Ls       = Ls;
    actionPar.mass     = 0.000678;
    actionPar.M5       = 1.8;
    actionPar.b        = 1.5;
    actionPar.c        = 0.5;
    actionPar.boundary = "1 1 1 1";
    actionPar.twist    = "0. 0. 0. 0.";
    application.createModule<MAction::MobiusDWF>("mdwf_l", actionPar);

    MAction::MobiusDWFPar actionParF = actionPar;
    actionParF.gauge = "gaugef";
    application.createModule<MAction::MobiusDWFF>("mdwf_lf", actionParF);

    // ------------------------------------------------------------------
    // Mixed-precision solver (registers both "mcg_l" and "mcg_l_subtract"),
    // matching production's mcg_l shape. outerGuesser is set below once the
    // eigenpack-backed guesser module exists -- see file header.
    // ------------------------------------------------------------------
    MSolver::MixedPrecisionRBPrecCGPar solverPar;
    solverPar.innerAction       = "mdwf_lf";
    solverPar.outerAction       = "mdwf_l";
    solverPar.maxInnerIteration = 2000;
    solverPar.maxOuterIteration = 2;
    solverPar.residual          = 1e-4;
    solverPar.innerGuesser      = "";
    solverPar.outerGuesser      = "guesser";
    solverPar.ifCGD             = false;
    application.createModule<MSolver::MixedPrecisionRBPrecCG>("mcg_l", solverPar);

    // ------------------------------------------------------------------
    // Time-diluted noise (same as noise_ud)
    // ------------------------------------------------------------------
    application.createModule<MNoise::TimeDilutedSpinColorDiagonal>("noise_ud");

    // ------------------------------------------------------------------
    // Synthetic eigenpack (stand-in for epackD)
    // ------------------------------------------------------------------
    MUtilities::RandomFermionEigenPackPar epackPar;
    epackPar.size    = Nl;
    epackPar.Ls      = Ls;
    epackPar.evalMin = evalMin;
    epackPar.evalMax = evalMax;
    application.createModule<MUtilities::RandomFermionEigenPack>("epackD", epackPar);

    // ------------------------------------------------------------------
    // Exact deflation guesser built from the same eigenpack, wired in as
    // mcg_l's outerGuesser. Production always configures a real deflation
    // guesser here; leaving it empty falls back to Grid's ZeroGuesser,
    // which (unlike DeflatedGuesser/ExactDeflatedGuesser) never sets
    // guess.Checkerboard() -- that mistagged zero field then survives into
    // SchurRedBlackBase's subGuess subtraction (active here since A2AVectors
    // auto-selects the "_subtract" solver whenever Nl_ > 0) and trips
    // Lattice_ET.h's checkerboard-consistency GRID_ASSERT. Wiring up the
    // guesser both mirrors production and avoids that unrelated latent bug.
    // ------------------------------------------------------------------
    MGuesser::ExactDeflationPar guesserPar;
    guesserPar.eigenPack = "epackD";
    guesserPar.size      = Nl;
    application.createModule<MGuesser::ExactDeflation>("guesser", guesserPar);

    // ------------------------------------------------------------------
    // A2A vectors (same module pattern as a2a_l)
    // ------------------------------------------------------------------
    MSolver::A2AVectorsPar a2aPar;
    a2aPar.noise     = "noise_ud";
    a2aPar.action    = "mdwf_l";
    a2aPar.eigenPack = "epackD";
    a2aPar.solver    = "mcg_l";
    a2aPar.output    = "";
    a2aPar.multiFile = true;
    application.createModule<MSolver::A2AVectors>("a2a_l", a2aPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
