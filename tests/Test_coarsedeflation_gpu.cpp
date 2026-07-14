// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_coarsedeflation_gpu.cpp
 *
 * Regression test for the 64I CGl memory redesign: checks that
 * MGuesser::CoarseDeflation200F (deflating directly against the compressed
 * coarse eigenpack, epack_coarse) reproduces the same guess -- and the same
 * downstream RBPrecCG solve -- as the current production path of
 * MUtilities::EigenPackLCDecompress200F + MGuesser::ExactDeflationF (dense
 * decompression into a resident single-precision array, epack).
 *
 * This test loads REAL 64I production data (the real coarse eigenpack and
 * the real gauge config it was computed from), unlike Test_a2a_evalfpe_gpu.cpp
 * which uses synthetic in-memory data. The filestem/file paths below are
 * placeholders -- fill in the real Frontier lustre paths before running.
 *
 * sizeFine is intentionally left at the full 200 (matches the local
 * coherence basis dimension baked into CoarseField's per-site type, see
 * Hadrons/EigenPack.hpp's CoarseFermionEigenPack alias) -- truncating it
 * would drop real coefficient contributions from every single reconstructed
 * vector, not just reduce the number of modes tested. sizeCoarse (the
 * number of independent approximate eigenmodes) is safely truncated to 50
 * for a fast smoke test; each mode reconstructs at full fidelity regardless
 * of how many others are loaded.
 *
 * A2AVectors is deliberately NOT used here -- only the guesser and the CG
 * solve are under test, so both RBPrecCGF instances are fed directly from a
 * single plain (undiluted) random source, generated once and shared between
 * both branches for a fair diff. Since neither the guesser nor RBPrecCGF
 * touch A2AVectors/SaveBinnedA2AVecs (the only reason production needs a
 * double-precision boundary at all), the entire test stays single precision
 * end to end -- no mixed-precision solver, no double eigenpack cast.
 *
 * Diff strategy, in order of how much each actually tells you:
 *   1. Diff the two guessers' raw output vectors directly (guess_diff),
 *      before any CG runs. This is the real regression signal: CG converges
 *      to the same fixed point regardless of initial guess quality (given
 *      enough iterations), so a converged final-solution diff alone would
 *      be a weak test of guesser correctness specifically.
 *   2. Iteration-count parity between the two RBPrecCGF solves -- since CG
 *      is deterministic, matching guesses should produce matching iteration
 *      counts. No special code needed: Grid's ConjugateGradient already
 *      logs this per solve.
 *   3. Diff the two final solution vectors (sol_diff) as a basic sanity
 *      check that both solves hit the same linear system correctly.
 *
 * Because deflation quality from only 50 (of 2000) coarse modes is weak,
 * both RBPrecCGF instances use a loose residual and a generous fixed
 * maxIteration as a backstop -- tune both down via --residual/--maxIteration
 * once guesser parity is confirmed.
 *
 * Usage (real 64I volume; loading the real coarse eigenpack forces the full
 * 64.64.64.128 volume regardless of node count, since that's the volume the
 * file was written at). Keep this to <=64 nodes; 32 nodes (256 ranks at 8
 * ranks/node) is enough since no A2A vectors are ever allocated:
 *   srun -N 32 ... ./Test_coarsedeflation_gpu --grid 64.64.64.128 \
 *       --mpi 4.4.4.4 --seed "1 2 3 4" --sizeCoarse 50 --maxIteration 1000 \
 *       --residual 1.0
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

// ---------------------------------------------------------------------------
// Test-local helper modules. None of these belong in production Hadrons --
// they exist only to drive/diff the guesser and solver objects directly,
// bypassing A2AVectors. Plain MODULE_REGISTER (not _TMP): these only exist
// in this test file, so there is no companion .cpp providing an explicit
// instantiation for the linker to find.
// ---------------------------------------------------------------------------
BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MUtilities)

// One random field on the fine red-black (Odd) grid -- the same grid the
// eigenpack's evec/guessers operate on. Used as the shared source for the
// raw guesser diff (guess_old/guess_new below).
class RandomOddFermionPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(RandomOddFermionPar,
                                    unsigned int, Ls);
};

template <typename FImpl>
class TRandomOddFermion: public Module<RandomOddFermionPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TRandomOddFermion(const std::string name)
    : Module<RandomOddFermionPar>(name)
    {}
    virtual ~TRandomOddFermion(void) {};
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
        envCreate(FermionField, getName(), par().Ls,
                  envGetRbGrid(FermionField, par().Ls));
    }
    virtual void execute(void)
    {
        // Same pattern as Test_a2a_evalfpe_gpu.cpp's TRandomFermionEigenPack:
        // GridParallelRNG::fill's automatic checkerboard handling needs a
        // full 5D grid to build its internal temporary against, so fill on
        // the full grid and project down rather than filling the RB field
        // directly.
        LOG(Message) << "Generating random field on the Odd checkerboard"
                     << std::endl;

        FermionField full(envGetGrid(FermionField, par().Ls));
        auto         &odd = envGet(FermionField, getName());

        random(rng4d(), full);
        pickCheckerboard(Odd, odd, full);
    }
};

MODULE_REGISTER(RandomOddFermionF, ARG(TRandomOddFermion<FIMPLF>), MUtilities);

// Calls a named LinearFunction<FermionField> (i.e. a Hadrons guesser module)
// directly on a named source field and stores the result -- no CG involved.
class CallGuesserPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CallGuesserPar,
                                    std::string,  guesser,
                                    std::string,  source,
                                    unsigned int, Ls);
};

template <typename FImpl>
class TCallGuesser: public Module<CallGuesserPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TCallGuesser(const std::string name)
    : Module<CallGuesserPar>(name)
    {}
    virtual ~TCallGuesser(void) {};
    virtual std::vector<std::string> getInput(void)
    {
        return {par().guesser, par().source};
    }
    virtual std::vector<std::string> getOutput(void)
    {
        return {getName()};
    }
    virtual void setup(void)
    {
        envCreate(FermionField, getName(), par().Ls,
                  envGetRbGrid(FermionField, par().Ls));
    }
    virtual void execute(void)
    {
        LOG(Message) << "Calling guesser '" << par().guesser
                     << "' on source '" << par().source << "'" << std::endl;

        auto &guesser = envGet(LinearFunction<FermionField>, par().guesser);
        auto &src     = envGet(FermionField, par().source);
        auto &guess   = envGet(FermionField, getName());

        guesser(src, guess);
    }
};

MODULE_REGISTER(CallGuesserF, ARG(TCallGuesser<FIMPLF>), MUtilities);

// Calls a named Solver<FermionField> (i.e. a Hadrons solver module, here
// RBPrecCGF) directly on a named full (non-checkerboarded) source field and
// stores the result.
class CallSolverPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CallSolverPar,
                                    std::string,  solver,
                                    std::string,  source,
                                    unsigned int, Ls);
};

template <typename FImpl>
class TCallSolver: public Module<CallSolverPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
public:
    TCallSolver(const std::string name)
    : Module<CallSolverPar>(name)
    {}
    virtual ~TCallSolver(void) {};
    virtual std::vector<std::string> getInput(void)
    {
        return {par().solver, par().source};
    }
    virtual std::vector<std::string> getOutput(void)
    {
        return {getName()};
    }
    virtual void setup(void)
    {
        envCreate(FermionField, getName(), par().Ls,
                  envGetGrid(FermionField, par().Ls));
    }
    virtual void execute(void)
    {
        LOG(Message) << "Calling solver '" << par().solver
                     << "' on source '" << par().source << "'" << std::endl;

        auto &solver = envGet(Solver, par().solver);
        auto &src    = envGet(FermionField, par().source);
        auto &sol    = envGet(FermionField, getName());

        solver(sol, src);
    }
};

MODULE_REGISTER(CallSolverF, ARG(TCallSolver<FIMPLF>), MUtilities);

// Logs |a|, |a - b| and |a - b|/|a| for two named fields of the same type.
// No product -- this is a terminal/sink module.
class FieldDiffPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(FieldDiffPar,
                                    std::string, a,
                                    std::string, b);
};

template <typename FImpl>
class TFieldDiff: public Module<FieldDiffPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TFieldDiff(const std::string name)
    : Module<FieldDiffPar>(name)
    {}
    virtual ~TFieldDiff(void) {};
    virtual std::vector<std::string> getInput(void)
    {
        return {par().a, par().b};
    }
    virtual std::vector<std::string> getOutput(void)
    {
        return {};
    }
    virtual void setup(void)
    {}
    virtual void execute(void)
    {
        auto &a = envGet(FermionField, par().a);
        auto &b = envGet(FermionField, par().b);

        FermionField diff(a.Grid());
        diff = a - b;

        RealD na    = std::sqrt(norm2(a));
        RealD nDiff = std::sqrt(norm2(diff));
        RealD rel   = (na > 0.) ? nDiff / na : nDiff;

        LOG(Message) << "FieldDiff '" << par().a << "' vs '" << par().b
                     << "':" << std::endl;
        LOG(Message) << "  |a|       = " << na    << std::endl;
        LOG(Message) << "  |a - b|   = " << nDiff << std::endl;
        LOG(Message) << "  |a-b|/|a| = " << rel   << std::endl;
    }
};

MODULE_REGISTER(FieldDiffF, ARG(TFieldDiff<FIMPLF>), MUtilities);

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
    globalPar.runId                    = "coarsedeflation_regression";
    // Modules below are created in dependency order, so the naive
    // (sequential, non-genetic) scheduler is valid here, same as
    // Test_a2a_evalfpe_gpu.cpp.
    globalPar.scheduler.schedulerType  = "naive";
    globalPar.genetic.maxGen           = 1000;
    globalPar.genetic.maxCstGen        = 200;
    globalPar.genetic.popSize          = 20;
    globalPar.genetic.mutationRate     = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments
    // ------------------------------------------------------------------
    unsigned int sizeCoarse   = 50;
    unsigned int maxIteration = 1000;
    double       residual     = 1.0;
    if (GridCmdOptionExists(argv, argv + argc, "--sizeCoarse"))
        sizeCoarse = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--sizeCoarse"));
    if (GridCmdOptionExists(argv, argv + argc, "--maxIteration"))
        maxIteration = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--maxIteration"));
    if (GridCmdOptionExists(argv, argv + argc, "--residual"))
        residual = std::stod(GridCmdOptionPayload(argv, argv + argc, "--residual"));

    const unsigned int Ls = 12;

    // ------------------------------------------------------------------
    // Real production gauge config (single-precision cast for the action
    // under test). Placeholder path -- fill in the real Frontier lustre
    // path before running.
    // ------------------------------------------------------------------
    MIO::LoadNerscPar gaugePar;
    gaugePar.file = "REPLACE_ME_gauge_config_path";
    application.createModule<MIO::LoadNersc>("gauge", gaugePar);

    MUtilities::PrecisionCastPar gaugefPar;
    gaugefPar.field = "gauge";
    application.createModule<MUtilities::GaugeSinglePrecisionCast>("gaugef", gaugefPar);

    // ------------------------------------------------------------------
    // Real production coarse eigenpack. sizeFine stays at the full 200
    // (see header comment); sizeCoarse is truncated for a fast test.
    // Placeholder filestem -- fill in the real Frontier lustre path before
    // running.
    // ------------------------------------------------------------------
    MIO::LoadCoarseEigenPackPar epackCoarsePar;
    epackCoarsePar.filestem      = "REPLACE_ME_coarse_eigenpack_filestem";
    epackCoarsePar.multiFile     = true;
    epackCoarsePar.sizeFine      = 200;
    epackCoarsePar.sizeCoarse    = sizeCoarse;
    epackCoarsePar.redBlack      = true;
    epackCoarsePar.Ls            = Ls;
    epackCoarsePar.blockSize     = "4 4 4 4 12";
    epackCoarsePar.orthogonalise = false;
    epackCoarsePar.gaugeXform    = "";
    application.createModule<MIO::LoadCoarseFermionEigenPack200F>("epack_coarse", epackCoarsePar);

    // ------------------------------------------------------------------
    // Old path: dense decompression into a resident single-precision array,
    // then exact deflation against that array -- exactly production's
    // current epack/guesser modules.
    // ------------------------------------------------------------------
    MUtilities::EigenPackLCDecompressPar epackPar;
    epackPar.epack      = "epack_coarse";
    epackPar.blockSize  = "4 4 4 4 12";
    epackPar.coarseSize = sizeCoarse;
    epackPar.Ls         = Ls;
    epackPar.output     = "";
    epackPar.multiFile  = true;
    application.createModule<MUtilities::EigenPackLCDecompress200F>("epack", epackPar);

    MGuesser::ExactDeflationPar guesserOldPar;
    guesserOldPar.eigenPack = "epack";
    guesserOldPar.size      = sizeCoarse;
    application.createModule<MGuesser::ExactDeflationF>("guesser_old", guesserOldPar);

    // ------------------------------------------------------------------
    // New path: deflate directly against the compressed coarse eigenpack,
    // no dense decompression at all.
    // ------------------------------------------------------------------
    MGuesser::CoarseDeflationPar guesserNewPar;
    guesserNewPar.eigenPack = "epack_coarse";
    guesserNewPar.size      = sizeCoarse;
    application.createModule<MGuesser::CoarseDeflation200F>("guesser_new", guesserNewPar);

    // ------------------------------------------------------------------
    // Single-precision Mobius DWF action, standard 64I parameters, matching
    // production's mdwff_l.
    // ------------------------------------------------------------------
    MAction::MobiusDWFPar actionParF;
    actionParF.gauge    = "gaugef";
    actionParF.Ls       = Ls;
    actionParF.mass     = 0.000678;
    actionParF.M5       = 1.8;
    actionParF.b        = 1.5;
    actionParF.c        = 0.5;
    actionParF.boundary = "1 1 1 -1";
    actionParF.twist    = "0. 0. 0. 0.";
    application.createModule<MAction::MobiusDWFF>("mdwff_l", actionParF);

    // ------------------------------------------------------------------
    // Shared random sources -- one per branch pair, generated once and
    // reused identically by both old/new modules for a fair diff.
    // ------------------------------------------------------------------
    MUtilities::RandomOddFermionPar srcOddPar;
    srcOddPar.Ls = Ls;
    application.createModule<MUtilities::RandomOddFermionF>("src_odd", srcOddPar);

    MUtilities::RandomFieldPar srcFullPar;
    srcFullPar.Ls = Ls;
    application.createModule<MUtilities::RandomFermionF>("src_full", srcFullPar);

    // ------------------------------------------------------------------
    // Primary diff: raw guesser output, no CG involved.
    // ------------------------------------------------------------------
    MUtilities::CallGuesserPar guessOldPar;
    guessOldPar.guesser = "guesser_old";
    guessOldPar.source  = "src_odd";
    guessOldPar.Ls      = Ls;
    application.createModule<MUtilities::CallGuesserF>("guess_old", guessOldPar);

    MUtilities::CallGuesserPar guessNewPar;
    guessNewPar.guesser = "guesser_new";
    guessNewPar.source  = "src_odd";
    guessNewPar.Ls      = Ls;
    application.createModule<MUtilities::CallGuesserF>("guess_new", guessNewPar);

    MUtilities::FieldDiffPar guessDiffPar;
    guessDiffPar.a = "guess_old";
    guessDiffPar.b = "guess_new";
    application.createModule<MUtilities::FieldDiffF>("guess_diff", guessDiffPar);

    // ------------------------------------------------------------------
    // Secondary/tertiary checks: full RBPrecCGF solves off the same
    // source, one per guesser. Loose residual + generous maxIteration
    // since 50-mode deflation converges slowly; both must actually
    // converge for the solution diff to mean anything.
    // ------------------------------------------------------------------
    MSolver::RBPrecCGPar cgOldPar;
    cgOldPar.action       = "mdwff_l";
    cgOldPar.maxIteration = maxIteration;
    cgOldPar.residual     = residual;
    cgOldPar.guesser      = "guesser_old";
    application.createModule<MSolver::RBPrecCGF>("cg_old", cgOldPar);

    MSolver::RBPrecCGPar cgNewPar;
    cgNewPar.action       = "mdwff_l";
    cgNewPar.maxIteration = maxIteration;
    cgNewPar.residual     = residual;
    cgNewPar.guesser      = "guesser_new";
    application.createModule<MSolver::RBPrecCGF>("cg_new", cgNewPar);

    MUtilities::CallSolverPar solOldPar;
    solOldPar.solver = "cg_old";
    solOldPar.source = "src_full";
    solOldPar.Ls     = Ls;
    application.createModule<MUtilities::CallSolverF>("sol_old", solOldPar);

    MUtilities::CallSolverPar solNewPar;
    solNewPar.solver = "cg_new";
    solNewPar.source = "src_full";
    solNewPar.Ls     = Ls;
    application.createModule<MUtilities::CallSolverF>("sol_new", solNewPar);

    MUtilities::FieldDiffPar solDiffPar;
    solDiffPar.a = "sol_old";
    solDiffPar.b = "sol_new";
    application.createModule<MUtilities::FieldDiffF>("sol_diff", solDiffPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
