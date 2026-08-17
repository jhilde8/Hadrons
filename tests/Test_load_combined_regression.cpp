/*
 * Test_load_combined_regression.cpp
 *
 * Regression test: does MIO::LoadCombinedA2AVecs correctly reproduce the
 * same array MIO::LoadBinnedA2AVecs would give when reading the exact same
 * underlying vector data, laid out the old way (one monolithic bin size
 * spanning low+high modes together, as vectors were stored before the
 * multi-hit rework) versus the new way (low modes and one hit's high
 * modes in separate bin sizes and files).
 *
 * Phase 1 (plain Grid, no Hadrons): generate w and v (656 random vectors
 * each, independent), slice each into low (400) + high (256), and write
 * each to disk twice -- once as the full 656-vector array in 4 bins of
 * 164 (old-style layout), once as low/high in their own files at bin
 * sizes 200/128 (new-style layout). The pack-and-write logic mirrors
 * MIO::SaveBinnedA2AVecs::execute() exactly, just called directly instead
 * of through a Module -- carving named sub-ranges out of a
 * Hadrons-resident array isn't something any existing module supports,
 * and there's no established pattern here for reaching into the
 * Environment mid-run to inject them, so this stays outside the Hadrons
 * graph entirely.
 *
 * Phase 2 (real Hadrons::Application): load w and v both ways
 * (LoadBinnedA2AVecs<164> old-style, LoadCombinedA2AVecs<200,128>
 * new-style), compute one A2AMesonField with both operands old-style and
 * one with both new-style, and write both out. Since the old-style files
 * contain the same 656 values in the same order as the new-style
 * concatenation, the two meson fields are directly comparable with the
 * unmodified diff_mf.py:
 *   python3 diff_mf.py mf_old_out.0 mf_combined_out.0
 *
 * Usage:
 *   mpirun -n 1 ./Test_load_combined_regression --grid 4.4.4.8 --mpi 1.1.1.1
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Hadrons/A2AVectors.hpp>

using namespace Grid;
using namespace Hadrons;

typedef FIMPL::FermionField FermionField;

template <int binSize>
void writeBinned(const std::string &filestem, std::vector<FermionField> &v, int traj)
{
    typedef iVector<typename FIMPL::SiteSpinor, binSize> SiteSpinorSet;

    assert(v.size() % binSize == 0);
    int Nb = (int)v.size() / binSize;
    std::vector<Lattice<SiteSpinorSet>> bvec(Nb, v[0].Grid());

    for (int i = 0; i < (int)v.size(); i += binSize)
        for (int j = 0; j < binSize; ++j)
            pokeLorentz(bvec[i / binSize], v[i + j], j);

    A2AVectorsIo::write(filestem, bvec, true, traj);
}

// Generates one role's (w or v) low+high split and writes both the
// old-style monolithic file and the new-style low/high files.
void makeRole(const std::string &tag, GridCartesian &grid, GridParallelRNG &pRNG,
              int nLow, int nHigh, int traj)
{
    int nTotal = nLow + nHigh;

    std::vector<FermionField> full(nTotal, &grid);
    for (auto &f : full) random(pRNG, f);

    std::vector<FermionField> low(nLow, &grid);
    std::vector<FermionField> high(nHigh, &grid);
    for (int i = 0; i < nLow;  i++) low[i]  = full[i];
    for (int i = 0; i < nHigh; i++) high[i] = full[nLow + i];

    writeBinned<164>("old_"  + tag, full, traj);
    writeBinned<200>("low_"  + tag, low,  traj);
    writeBinned<128>("high_" + tag, high, traj);
}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    const int traj  = 0;
    const int nLow   = 400;
    const int nHigh  = 256;
    const int nTotal = nLow + nHigh;   // 656

    // ------------------------------------------------------------------
    // Phase 1: generate + slice + write, plain Grid, no Hadrons.
    // ------------------------------------------------------------------
    Coordinate latt_size   = GridDefaultLatt();
    Coordinate simd_layout = GridDefaultSimd(4, vComplexD::Nsimd());
    Coordinate mpi_layout  = GridDefaultMpi();
    GridCartesian grid(latt_size, simd_layout, mpi_layout);

    GridParallelRNG pRNG(&grid);
    pRNG.SeedFixedIntegers({1, 2, 3, 4});

    makeRole("w", grid, pRNG, nLow, nHigh, traj);
    makeRole("v", grid, pRNG, nLow, nHigh, traj);

    // ------------------------------------------------------------------
    // Phase 2: real Hadrons application -- load w and v both ways,
    // compute meson fields for comparison with diff_mf.py.
    // ------------------------------------------------------------------
    Application application;

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = traj;
    globalPar.trajCounter.end      = traj + 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "load_combined_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    auto addOld = [&](const std::string &tag)
    {
        MIO::LoadBinnedA2AVecsPar par;
        par.filestem  = "old_" + tag;
        par.multiFile = true;
        par.size      = nTotal;
        application.createModule<MIO::TLoadBinnedA2AVecs<FIMPL, 164>>(tag + "_old", par);
    };

    auto addCombined = [&](const std::string &tag)
    {
        MIO::LoadCombinedA2AVecsPar par;
        par.lowFilestem    = "low_" + tag;
        par.nLow           = nLow;
        par.highStem       = "high_" + tag;
        par.highExtensions = {""};
        par.nHighEach      = nHigh;
        application.createModule<MIO::TLoadCombinedA2AVecs<FIMPL, 200, 128>>(tag + "_combined", par);
    };

    addOld("w");
    addOld("v");
    addCombined("w");
    addCombined("v");

    MContraction::A2AMesonFieldPar mfOldPar;
    mfOldPar.block      = 128;
    mfOldPar.cacheBlock = 16;
    mfOldPar.left       = "w_old";
    mfOldPar.right      = "v_old";
    mfOldPar.output     = "mf_old_out";
    mfOldPar.gammas     = "Gamma5";
    mfOldPar.mom        = {"0 0 0"};
    application.createModule<MContraction::A2AMesonField>("mf_old", mfOldPar);

    MContraction::A2AMesonFieldPar mfCombPar = mfOldPar;
    mfCombPar.left   = "w_combined";
    mfCombPar.right  = "v_combined";
    mfCombPar.output = "mf_combined_out";
    application.createModule<MContraction::A2AMesonField>("mf_combined", mfCombPar);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
