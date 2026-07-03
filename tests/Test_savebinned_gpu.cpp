// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_savebinned_gpu.cpp
 *
 * Hadrons application: generates random A2A vectors and writes them to disk
 * with MIO::SaveBinnedA2AVecs (GPU path, pack via pokeLorentz) to check that
 * the module runs and produces the expected SciDAC output in the GPU build.
 *
 * binSize is fixed to 197 to match the binning used for the light-quark low
 * modes in production (see par.CGl.905.xml, module SaveBinnedA2AVecs197).
 * Swap the MIO::SaveBinnedA2AVecs197 type below for another registered size
 * (48, 64, 96, 128, 144, 149, 173, 192, 196, 216, 221, 298) to test a
 * different bin.
 *
 * Usage:
 *   mpirun -n 1 ./Test_savebinned_gpu --grid 4.4.4.8 --mpi 1.1.1.1 --seed "1 2 3 4"
 *
 * Output files: savebinned_gpu_out.0/elem0.bin, elem1.bin (multiFile=true).
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
    // Global parameters
    // ------------------------------------------------------------------
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = 0;
    globalPar.trajCounter.end      = 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "savebinned_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    // ------------------------------------------------------------------
    // Parse optional CLI arguments
    // ------------------------------------------------------------------
    const unsigned int binSize = 197;
    unsigned int       Nb      = 2;
    if (GridCmdOptionExists(argv, argv + argc, "--Nb"))
        Nb = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--Nb"));

    // ------------------------------------------------------------------
    // Random A2A vectors: Nb bins of binSize vectors each.
    // ------------------------------------------------------------------
    MUtilities::RandomVectorsPar rv;
    rv.size      = Nb * binSize;
    rv.Ls        = 1;
    rv.output    = "";
    rv.multiFile = false;

    application.createModule<MUtilities::RandomFermions>("a2a_test", rv);

    // ------------------------------------------------------------------
    // SaveBinnedA2AVecs - GPU path (pokeLorentz pack)
    // ------------------------------------------------------------------
    MIO::SaveBinnedA2AVecsPar savePar;
    savePar.filestem  = "savebinned_gpu_out";
    savePar.multiFile = true;
    savePar.a2a       = "a2a_test";

    application.createModule<MIO::SaveBinnedA2AVecs197>("save_test", savePar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
