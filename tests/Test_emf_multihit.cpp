/*
 * Test_emf_multihit.cpp
 *
 * Regression test: does A2AExtendedMesonField, fed a multi-hit combined-
 * format A2A vector array for its quark loop (MIO::LoadCombinedA2AVecs +
 * MUtilities::VectorPackRefSlice, pointer-referencing a shared low block
 * plus a specific hit's high block out of one resident array), reproduce
 * the same result as the old way -- a single dedicated per-hit file,
 * loaded and pointer-referenced in its entirety (the "reference to the
 * environment object" production used before this rework).
 *
 * Two things are exercised simultaneously:
 *
 *   1. External legs (left = l_v, right = s_v) carry the same hit
 *      degree of freedom as the MF regression test: old-style per-hit
 *      files (light: low duplicated in; strange: nothing to duplicate)
 *      vs. one new-style combined array spanning all hits. The combined
 *      EMF's N_i x N_j output block-decomposes into external hit-pair
 *      sub-blocks, each comparable to a separately-computed old-style
 *      pairwise EMF.
 *
 *   2. The quark loop (loop_vw1/loop_vw2) is fixed at a single hit index
 *      -loopHit- for the whole run, matching how production actually
 *      runs (one job per loop hit, "indexed by loop index"). The new
 *      path builds loop pointers via VectorPackRefSlice with segments
 *      [0,Nl) + [Nl+loopHit*Nh, Nl+(loopHit+1)*Nh) out of the *same*
 *      combined array used for the external leg -- non-contiguous, and
 *      only trivially equivalent to a full reference when loopHit = 0.
 *      The old path builds loop pointers via VectorPackRefSlice with one
 *      segment spanning the entirety of a dedicated old-style hit-
 *      -loopHit- file. Both feed the same, unmodified
 *      A2AExtendedMesonField module, which always calls LoopPropagatorPtr
 *      (bare LoopPropagator vs LoopPropagatorPtr equivalence is already
 *      covered at the Grid level by Test_extended_meson_field.cc).
 *
 * Two loop flavors are covered: light (loop_vw1/2 = l_v/l_w, low+high)
 * and strange (loop_vw1/2 = s_v/s_w, high-only). l_v/s_v double as both
 * the external legs and the loop's V leg -- this reuse is the entire
 * point of building the loop from pointers into a resident array rather
 * than a dedicated load. l_w/s_w are generated purely for the loop's W
 * leg, with the same dual old/new structure so the same offset
 * reconstruction gets exercised on both loop legs.
 *
 * Usage:
 *   mpirun -n 1 ./Test_emf_multihit --grid 4.4.4.8 --mpi 1.1.1.1
 *   mpirun -n 1 ./Test_emf_multihit --grid 4.4.4.8 --mpi 1.1.1.1 --skip-gen
 *   mpirun -n 1 ./Test_emf_multihit --grid 4.4.4.8 --mpi 1.1.1.1 --loop-hit 1
 *   mpirun -n 1 ./Test_emf_multihit --grid 4.4.4.8 --mpi 1.1.1.1 --light-hits 1 --strange-hits 3
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
    const int Nl   = 100;   // shared low modes, light only
    const int Nh   = 64;    // high modes per hit, both flavors

    int nLightHits   = 2;
    int nStrangeHits = 2;
    int loopHit      = 0;
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
    if (GridCmdOptionExists(argv, argv + argc, "--loop-hit"))
    {
        std::string arg = GridCmdOptionPayload(argv, argv + argc, "--loop-hit");
        GridCmdOptionInt(arg, loopHit);
    }
    bool skipGen = GridCmdOptionExists(argv, argv + argc, "--skip-gen");

    assert(loopHit >= 0 && loopHit < nLightHits && loopHit < nStrangeHits);

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

        // l_v: shared low modes, per-hit high modes, old-style files with
        // the low block duplicated in.
        auto lvLow = genRandom(Nl, grid, pRNG);
        writeBinned<Nl>("low_l_v", lvLow, traj);
        for (int k = 0; k < nLightHits; ++k)
        {
            std::string ext = std::to_string(k);
            auto high = genRandom(Nh, grid, pRNG);
            writeBinned<Nh>("high_l_v_hit" + ext, high, traj);

            std::vector<FermionField> old(Nl + Nh, &grid);
            for (int i = 0; i < Nl; ++i) old[i]      = lvLow[i];
            for (int i = 0; i < Nh; ++i) old[Nl + i] = high[i];
            writeBinned<Nl + Nh>("l_v_old_hit" + ext, old, traj);
        }

        // l_w: same structure, independently random -- needed purely for
        // the loop's W leg.
        auto lwLow = genRandom(Nl, grid, pRNG);
        writeBinned<Nl>("low_l_w", lwLow, traj);
        for (int k = 0; k < nLightHits; ++k)
        {
            std::string ext = std::to_string(k);
            auto high = genRandom(Nh, grid, pRNG);
            writeBinned<Nh>("high_l_w_hit" + ext, high, traj);

            std::vector<FermionField> old(Nl + Nh, &grid);
            for (int i = 0; i < Nl; ++i) old[i]      = lwLow[i];
            for (int i = 0; i < Nh; ++i) old[Nl + i] = high[i];
            writeBinned<Nl + Nh>("l_w_old_hit" + ext, old, traj);
        }

        // s_v: no low block -- the high file doubles as the old-style
        // per-hit file.
        for (int k = 0; k < nStrangeHits; ++k)
        {
            auto high = genRandom(Nh, grid, pRNG);
            writeBinned<Nh>("s_v_hit" + std::to_string(k), high, traj);
        }

        // s_w: same, independently random.
        for (int k = 0; k < nStrangeHits; ++k)
        {
            auto high = genRandom(Nh, grid, pRNG);
            writeBinned<Nh>("s_w_hit" + std::to_string(k), high, traj);
        }
    }
    else
    {
        LOG(Message) << "--skip-gen: reusing existing files on disk" << std::endl;
    }

    // ------------------------------------------------------------------
    // Phase 2: real Hadrons application.
    // ------------------------------------------------------------------
    Application application;

    Application::GlobalPar globalPar;
    globalPar.trajCounter.start    = traj;
    globalPar.trajCounter.end      = traj + 1;
    globalPar.trajCounter.step     = 1;
    globalPar.runId                = "emf_multihit_regression";
    globalPar.genetic.maxGen       = 1000;
    globalPar.genetic.maxCstGen    = 200;
    globalPar.genetic.popSize      = 20;
    globalPar.genetic.mutationRate = .1;
    application.setPar(globalPar);

    auto addLightCombined = [&](const std::string &role)
    {
        MIO::LoadCombinedA2AVecsPar par;
        par.lowFilestem = "low_l_" + role;
        par.nLow        = Nl;
        par.highStem    = "high_l_" + role + "_hit";
        for (int k = 0; k < nLightHits; ++k) par.highExtensions.push_back(std::to_string(k));
        par.nHighEach   = Nh;

        std::string name = "l_" + role + "_combined";
        application.createModule<MIO::TLoadCombinedA2AVecs<FIMPL, Nl, Nh>>(name, par);

        return name;
    };

    auto addStrangeCombined = [&](const std::string &role)
    {
        MIO::LoadCombinedA2AVecsPar par;
        par.nLow      = 0;
        par.highStem  = "s_" + role + "_hit";
        for (int k = 0; k < nStrangeHits; ++k) par.highExtensions.push_back(std::to_string(k));
        par.nHighEach = Nh;

        std::string name = "s_" + role + "_combined";
        application.createModule<MIO::TLoadCombinedA2AVecs<FIMPL, Nl, Nh>>(name, par);

        return name;
    };

    auto addLightOld = [&](const std::string &role, int k)
    {
        MIO::LoadBinnedA2AVecsPar par;
        par.filestem  = "l_" + role + "_old_hit" + std::to_string(k);
        par.multiFile = true;
        par.size      = Nl + Nh;

        std::string name = "l_" + role + "_old_hit" + std::to_string(k);
        application.createModule<MIO::TLoadBinnedA2AVecs<FIMPL, Nl + Nh>>(name, par);

        return name;
    };

    auto addStrangeOld = [&](const std::string &role, int k)
    {
        MIO::LoadBinnedA2AVecsPar par;
        par.filestem  = "s_" + role + "_hit" + std::to_string(k);
        par.multiFile = true;
        par.size      = Nh;

        std::string name = "s_" + role + "_old_hit" + std::to_string(k);
        application.createModule<MIO::TLoadBinnedA2AVecs<FIMPL, Nh>>(name, par);

        return name;
    };

    std::string lightVComb   = addLightCombined("v");
    std::string lightWComb   = addLightCombined("w");
    std::string strangeVComb = addStrangeCombined("v");
    std::string strangeWComb = addStrangeCombined("w");

    std::vector<std::string> lightVOld(nLightHits), lightWOld(nLightHits);
    for (int k = 0; k < nLightHits; ++k)
    {
        lightVOld[k] = addLightOld("v", k);
        lightWOld[k] = addLightOld("w", k);
    }

    std::vector<std::string> strangeVOld(nStrangeHits), strangeWOld(nStrangeHits);
    for (int k = 0; k < nStrangeHits; ++k)
    {
        strangeVOld[k] = addStrangeOld("v", k);
        strangeWOld[k] = addStrangeOld("w", k);
    }

    // Loop pointer construction, fixed at hit index loopHit for the
    // whole run. New: offset into the combined array (non-contiguous
    // for loopHit > 0). Old: one trivial segment spanning the whole
    // dedicated hit-loopHit file.
    auto addPackSlice = [&](const std::string &name, const std::string &source,
                             const std::vector<unsigned int> &offsets,
                             const std::vector<unsigned int> &counts)
    {
        MUtilities::VectorPackRefSlicePar par;
        par.source  = source;
        par.offsets = offsets;
        par.counts  = counts;
        application.createModule<MUtilities::TVectorPackRefSlice<FIMPL::FermionField>>(name, par);
    };

    std::string kStr = std::to_string(loopHit);
    unsigned int uNl = (unsigned int)Nl;
    unsigned int uNh = (unsigned int)Nh;
    unsigned int highOffset = (unsigned int)(Nl + loopHit * Nh);
    unsigned int strangeOffset = (unsigned int)(loopHit * Nh);

    std::string loop1LightNew = "loop_l_v_new_hit" + kStr;
    addPackSlice(loop1LightNew, lightVComb, {0u, highOffset}, {uNl, uNh});
    std::string loop2LightNew = "loop_l_w_new_hit" + kStr;
    addPackSlice(loop2LightNew, lightWComb, {0u, highOffset}, {uNl, uNh});

    std::string loop1LightOld = "loop_l_v_old_hit" + kStr;
    addPackSlice(loop1LightOld, lightVOld[loopHit], {0u}, {uNl + uNh});
    std::string loop2LightOld = "loop_l_w_old_hit" + kStr;
    addPackSlice(loop2LightOld, lightWOld[loopHit], {0u}, {uNl + uNh});

    std::string loop1StrangeNew = "loop_s_v_new_hit" + kStr;
    addPackSlice(loop1StrangeNew, strangeVComb, {strangeOffset}, {uNh});
    std::string loop2StrangeNew = "loop_s_w_new_hit" + kStr;
    addPackSlice(loop2StrangeNew, strangeWComb, {strangeOffset}, {uNh});

    std::string loop1StrangeOld = "loop_s_v_old_hit" + kStr;
    addPackSlice(loop1StrangeOld, strangeVOld[loopHit], {0u}, {uNh});
    std::string loop2StrangeOld = "loop_s_w_old_hit" + kStr;
    addPackSlice(loop2StrangeOld, strangeWOld[loopHit], {0u}, {uNh});

    auto addEmf = [&](const std::string &name, const std::string &left, const std::string &right,
                       const std::string &loop1, const std::string &loop2)
    {
        MContraction::A2AExtendedMesonFieldPar par;
        par.block      = 64;
        par.cacheBlock = 16;
        par.types      = "0";
        par.left       = left;
        par.right      = right;
        par.loop_vw1   = loop1;
        par.loop_vw2   = loop2;
        par.output     = name + "_out";
        par.gammas1    = "Identity";
        par.gammas2    = "Identity";
        application.createModule<MContraction::TA2AExtendedMesonField<FIMPL>>(name, par);
    };

    // Light-loop: one combined-array EMF, plus one old-style EMF per
    // external hit pair, all sharing the same fixed-at-loopHit loop.
    addEmf("mf_light_loop_combined", lightVComb, strangeVComb, loop1LightNew, loop2LightNew);
    for (int i = 0; i < nLightHits; ++i)
        for (int j = 0; j < nStrangeHits; ++j)
            addEmf("mf_light_loop_old_" + std::to_string(i) + "_" + std::to_string(j),
                   lightVOld[i], strangeVOld[j], loop1LightOld, loop2LightOld);

    // Strange-loop: same pattern.
    addEmf("mf_strange_loop_combined", lightVComb, strangeVComb, loop1StrangeNew, loop2StrangeNew);
    for (int i = 0; i < nLightHits; ++i)
        for (int j = 0; j < nStrangeHits; ++j)
            addEmf("mf_strange_loop_old_" + std::to_string(i) + "_" + std::to_string(j),
                   lightVOld[i], strangeVOld[j], loop1StrangeOld, loop2StrangeOld);

    application.run();

    Grid_finalize();
    return EXIT_SUCCESS;
}
