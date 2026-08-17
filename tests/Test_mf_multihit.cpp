/*
 * Test_mf_multihit.cpp
 *
 * Regression test: does A2AMesonField, fed multi-hit combined-format A2A
 * vectors (MIO::LoadCombinedA2AVecs, pointer-referencing a shared low
 * block + per-hit high blocks), reproduce the same result as computing
 * it the old way -- one full (low+high) vector array per hit, with the
 * low modes literally duplicated into every hit's file, loaded via
 * MIO::LoadBinnedA2AVecs and contracted pairwise, hit by hit.
 *
 * Covers all four flavor combinations that occur in this project's
 * meson fields: light-light (low modes on both operands), light-strange
 * and strange-light (low modes on one operand only), and strange-strange
 * (no low modes on either operand). Light and strange hit counts are
 * independent CLI parameters, so the 1-hit collapse can be exercised on
 * either side without recompiling.
 *
 * Phase 1 (plain Grid, no Hadrons): generate nLightHits/nStrangeHits
 * independent sets of random high-mode vectors per role (light_w,
 * light_v, strange_w, strange_v), plus one shared set of low-mode
 * vectors per light role. Write:
 *   - new-style files: one low file per light role, one high file per
 *     hit per role (light and strange alike).
 *   - old-style files: for light, one monolithic per-hit file per role
 *     with the shared low vectors duplicated in ahead of that hit's high
 *     vectors; for strange, the high file itself doubles as the
 *     old-style file (nothing to duplicate -- no low block exists).
 *
 * Phase 2 (real Hadrons::Application): load combined-format light/
 * strange w/v (one call per role, spanning all hits), and separately
 * load old-style w/v per hit. Compute one A2AMesonField per flavor
 * combination from the combined-format operands (mf_{TYPE}_combined_out,
 * the full block matrix in one shot), and one A2AMesonField per hit pair
 * per flavor combination from the old-style operands
 * (mf_{TYPE}_old_{i}_{j}_out, the reference sub-blocks). A future diff
 * script slices the combined-format matrices into blocks and compares
 * them against the corresponding old-style outputs.
 *
 * Usage:
 *   mpirun -n 1 ./Test_mf_multihit --grid 4.4.4.8 --mpi 1.1.1.1
 *   mpirun -n 1 ./Test_mf_multihit --grid 4.4.4.8 --mpi 1.1.1.1 --skip-gen
 *   mpirun -n 1 ./Test_mf_multihit --grid 4.4.4.8 --mpi 1.1.1.1 --light-hits 1 --strange-hits 3
 */

#include <Hadrons/Global.hpp>
#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include <Hadrons/A2AVectors.hpp>

using namespace Grid;
using namespace Hadrons;

typedef FIMPL::FermionField FermionField;

template <int binSize>
void writeBinned(const std::string &filestem, std::vector<FermionField> &v, int traj)
{
    typedef typename FIMPL::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;

    assert(v.size() % binSize == 0);
    int Nb = (int)v.size() / binSize;
    std::vector<Lattice<SiteSpinorSet>> bvec(Nb, v[0].Grid());

    for (int i = 0; i < (int)v.size(); i += binSize)
        for (int j = 0; j < binSize; ++j)
            pokeLorentz(bvec[i / binSize], v[i + j], j);

    A2AVectorsIo::write(filestem, bvec, true, traj);
}

std::vector<FermionField> genRandom(int n, GridCartesian &grid, GridParallelRNG &pRNG)
{
    std::vector<FermionField> v(n, &grid);
    for (auto &f : v) random(pRNG, f);

    return v;
}

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    const int traj = 0;
    const int Nl   = 100;   // light low modes, shared across hits
    const int Nh   = 64;    // high modes per hit, both flavors

    int nLightHits   = 2;
    int nStrangeHits = 2;
    if (GridCmdOptionExists(argv, argv + argc, "--light-hits"))
    {
        std::string arg = GridCmdOptionPayload(argv, argv + argc, "--light-hits");
        GridCmdOptionInt(arg, nLightHits);
    }
    if (GridCmdOptionExists(argv, argv + argc, "--strange-hits"))
    {
        std::string arg = GridCmdOptionPayload(argv, argv + argc, "--strange-hits");
        GridCmdOptionInt(arg, nStrangeHits);
    }
    bool skipGen = GridCmdOptionExists(argv, argv + argc, "--skip-gen");

    // ------------------------------------------------------------------
    // Phase 1: generate + write, plain Grid, no Hadrons. Skipped with
    // --skip-gen, reusing whatever files are already on disk.
    // ------------------------------------------------------------------
    if (!skipGen)
    {
        Coordinate latt_size   = GridDefaultLatt();
        Coordinate simd_layout = GridDefaultSimd(4, vComplexD::Nsimd());
        Coordinate mpi_layout  = GridDefaultMpi();
        GridCartesian grid(latt_size, simd_layout, mpi_layout);

        GridParallelRNG pRNG(&grid);
        pRNG.SeedFixedIntegers({1, 2, 3, 4});

        // Shared light low modes, new-style files.
        auto lightLow_w = genRandom(Nl, grid, pRNG);
        auto lightLow_v = genRandom(Nl, grid, pRNG);
        writeBinned<Nl>("low_light_w", lightLow_w, traj);
        writeBinned<Nl>("low_light_v", lightLow_v, traj);

        // Light: per-hit high modes (new-style), plus old-style files
        // with the shared low block duplicated in ahead of the high part.
        for (int k = 0; k < nLightHits; ++k)
        {
            std::string ext = std::to_string(k);
            auto highW = genRandom(Nh, grid, pRNG);
            auto highV = genRandom(Nh, grid, pRNG);

            writeBinned<Nh>("high_light_w_hit" + ext, highW, traj);
            writeBinned<Nh>("high_light_v_hit" + ext, highV, traj);

            std::vector<FermionField> oldW(Nl + Nh, &grid);
            for (int i = 0; i < Nl; ++i) oldW[i]      = lightLow_w[i];
            for (int i = 0; i < Nh; ++i) oldW[Nl + i] = highW[i];
            writeBinned<Nl + Nh>("old_light_w_hit" + ext, oldW, traj);

            std::vector<FermionField> oldV(Nl + Nh, &grid);
            for (int i = 0; i < Nl; ++i) oldV[i]      = lightLow_v[i];
            for (int i = 0; i < Nh; ++i) oldV[Nl + i] = highV[i];
            writeBinned<Nl + Nh>("old_light_v_hit" + ext, oldV, traj);
        }

        // Strange: no low block on either side, so the high file itself
        // is already in old-style form -- one write per hit, used both
        // ways in Phase 2.
        for (int k = 0; k < nStrangeHits; ++k)
        {
            std::string ext = std::to_string(k);
            auto highW = genRandom(Nh, grid, pRNG);
            auto highV = genRandom(Nh, grid, pRNG);

            writeBinned<Nh>("strange_w_hit" + ext, highW, traj);
            writeBinned<Nh>("strange_v_hit" + ext, highV, traj);
        }
    }
    else
    {
        LOG(Message) << "--skip-gen: reusing existing files on disk" << std::endl;
    }

    // ------------------------------------------------------------------
    // Phase 2: real Hadrons application -- load both ways, compute
    // combined-format and old-style-pairwise meson fields for every
    // flavor combination.
    // ------------------------------------------------------------------
    Application application;

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = traj;
    globalPar.trajCounter.end      = traj + 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "mf_multihit_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    auto addLightCombined = [&](const std::string &tag)
    {
        MIO::LoadCombinedA2AVecsPar par;
        par.lowFilestem = "low_light_" + tag;
        par.nLow        = Nl;
        par.highStem    = "high_light_" + tag + "_hit";
        for (int k = 0; k < nLightHits; ++k) par.highExtensions.push_back(std::to_string(k));
        par.nHighEach   = Nh;

        std::string name = "light_" + tag + "_combined";
        application.createModule<MIO::TLoadCombinedA2AVecs<FIMPL, Nl, Nh>>(name, par);

        return name;
    };

    auto addStrangeCombined = [&](const std::string &tag)
    {
        MIO::LoadCombinedA2AVecsPar par;
        par.nLow      = 0;
        par.highStem  = "strange_" + tag + "_hit";
        for (int k = 0; k < nStrangeHits; ++k) par.highExtensions.push_back(std::to_string(k));
        par.nHighEach = Nh;

        std::string name = "strange_" + tag + "_combined";
        application.createModule<MIO::TLoadCombinedA2AVecs<FIMPL, Nl, Nh>>(name, par);

        return name;
    };

    auto addLightOld = [&](const std::string &tag, int k)
    {
        MIO::LoadBinnedA2AVecsPar par;
        par.filestem  = "old_light_" + tag + "_hit" + std::to_string(k);
        par.multiFile = true;
        par.size      = Nl + Nh;

        std::string name = "light_" + tag + "_old_" + std::to_string(k);
        application.createModule<MIO::TLoadBinnedA2AVecs<FIMPL, Nl + Nh>>(name, par);

        return name;
    };

    auto addStrangeOld = [&](const std::string &tag, int k)
    {
        MIO::LoadBinnedA2AVecsPar par;
        par.filestem  = "strange_" + tag + "_hit" + std::to_string(k);
        par.multiFile = true;
        par.size      = Nh;

        std::string name = "strange_" + tag + "_old_" + std::to_string(k);
        application.createModule<MIO::TLoadBinnedA2AVecs<FIMPL, Nh>>(name, par);

        return name;
    };

    std::string lightWComb   = addLightCombined("w");
    std::string lightVComb   = addLightCombined("v");
    std::string strangeWComb = addStrangeCombined("w");
    std::string strangeVComb = addStrangeCombined("v");

    std::vector<std::string> lightWOld(nLightHits), lightVOld(nLightHits);
    for (int k = 0; k < nLightHits; ++k)
    {
        lightWOld[k] = addLightOld("w", k);
        lightVOld[k] = addLightOld("v", k);
    }

    std::vector<std::string> strangeWOld(nStrangeHits), strangeVOld(nStrangeHits);
    for (int k = 0; k < nStrangeHits; ++k)
    {
        strangeWOld[k] = addStrangeOld("w", k);
        strangeVOld[k] = addStrangeOld("v", k);
    }

    auto addMf = [&](const std::string &name, const std::string &left, const std::string &right)
    {
        MContraction::A2AMesonFieldPar par;
        par.block      = 64;
        par.cacheBlock = 16;
        par.left       = left;
        par.right      = right;
        par.output     = name + "_out";
        par.gammas     = "Gamma5";
        par.mom        = {"0 0 0"};
        application.createModule<MContraction::A2AMesonField>(name, par);
    };

    // New-style combined meson fields -- one per flavor combination, the
    // full block matrix in one shot.
    addMf("mf_LL_combined", lightWComb,   lightVComb);
    addMf("mf_LS_combined", lightWComb,   strangeVComb);
    addMf("mf_SL_combined", strangeWComb, lightVComb);
    addMf("mf_SS_combined", strangeWComb, strangeVComb);

    // Old-style pairwise reference meson fields -- one per hit pair per
    // flavor combination, supplying the reference sub-blocks.
    for (int i = 0; i < nLightHits; ++i)
        for (int j = 0; j < nLightHits; ++j)
            addMf("mf_LL_old_" + std::to_string(i) + "_" + std::to_string(j),
                  lightWOld[i], lightVOld[j]);

    for (int i = 0; i < nLightHits; ++i)
        for (int j = 0; j < nStrangeHits; ++j)
            addMf("mf_LS_old_" + std::to_string(i) + "_" + std::to_string(j),
                  lightWOld[i], strangeVOld[j]);

    for (int i = 0; i < nStrangeHits; ++i)
        for (int j = 0; j < nLightHits; ++j)
            addMf("mf_SL_old_" + std::to_string(i) + "_" + std::to_string(j),
                  strangeWOld[i], lightVOld[j]);

    for (int i = 0; i < nStrangeHits; ++i)
        for (int j = 0; j < nStrangeHits; ++j)
            addMf("mf_SS_old_" + std::to_string(i) + "_" + std::to_string(j),
                  strangeWOld[i], strangeVOld[j]);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
