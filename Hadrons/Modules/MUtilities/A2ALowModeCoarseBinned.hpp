/*
 * A2ALowModeCoarseBinned.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MUtilities_A2ALowModeCoarseBinned_hpp_
#define Hadrons_MUtilities_A2ALowModeCoarseBinned_hpp_

#include <memory>

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Stream low-mode A2A V & W vectors straight to disk, binned, directly     *
 *  from a compressed (LC) coarse eigenpack. Low modes are decompressed/     *
 *  promoted on the fly, one bin at a time, so the full Nl-sized V/W arrays  *
 *  are never resident -- only one bin's worth of fields (compile-time       *
 *  binSize) at once. Coarse-deflation companion to                          *
 *  MUtilities::A2ALowModeBinned, which does the same streaming from an      *
 *  uncompressed (exact) eigenpack. Lives in MUtilities rather than MSolver  *
 *  because no solve ever happens here: low modes are built algebraically    *
 *  from the eigenvectors, so unlike MSolver::A2AHighModeVBinned there is    *
 *  no solver dependency at all.                                             *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ALowModeCoarseBinnedPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALowModeCoarseBinnedPar,
                                    std::string, eigenPack,
                                    std::string, action,
                                    std::string, output,
                                    std::string, schurConvention,
                                    unsigned int, checkInterval);
};

template <typename FImpl, typename FImplPack, int nBasis, int binSize>
class TA2ALowModeCoarseBinned : public Module<A2ALowModeCoarseBinnedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    FERM_TYPE_ALIASES(FImplPack, Pack);
    typedef CoarseFermionEigenPack<FImplPack, nBasis> EPack;
    typedef typename FImpl::SiteSpinor::vector_type   vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2ALowModeCoarseBinned(const std::string name);
    // destructor
    virtual ~TA2ALowModeCoarseBinned(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
    unsigned int Nb_{0};
    // Chosen from par().schurConvention at runtime in setup() -- see
    // MSolver::MixedPrecisionRBPrecCG's SOLVER_BODY macro and
    // Hadrons::A2AVectorsSchurBase for why this has to be a pointer to the
    // shared base class rather than a concrete stack object. Also used
    // directly (via a2a_->op()) for the in-program eigenvector check below,
    // so the check is guaranteed to validate evec_i against the exact same
    // operator used to build V/W, not a second, independently-configured one.
    std::unique_ptr<A2AVectorsSchurBase<FImpl>> a2a_;
};

MODULE_REGISTER_TMP(A2ALowModeCoarseBinned200Bin200,
    ARG(TA2ALowModeCoarseBinned<FIMPL, FIMPLF, 200, 200>), MUtilities);
MODULE_REGISTER_TMP(A2ALowModeCoarseBinned200Bin100,
    ARG(TA2ALowModeCoarseBinned<FIMPL, FIMPLF, 200, 100>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModeCoarseBinned200Bin200,
    ARG(TA2ALowModeCoarseBinned<ZFIMPL, ZFIMPLF, 200, 200>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModeCoarseBinned200Bin100,
    ARG(TA2ALowModeCoarseBinned<ZFIMPL, ZFIMPLF, 200, 100>), MUtilities);

/******************************************************************************
 *                   TA2ALowModeCoarseBinned implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
TA2ALowModeCoarseBinned<FImpl, FImplPack, nBasis, binSize>::TA2ALowModeCoarseBinned(const std::string name)
: Module<A2ALowModeCoarseBinnedPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
std::vector<std::string> TA2ALowModeCoarseBinned<FImpl, FImplPack, nBasis, binSize>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack, par().action};

    return in;
}

template <typename FImpl, typename FImplPack, int nBasis, int binSize>
std::vector<std::string> TA2ALowModeCoarseBinned<FImpl, FImplPack, nBasis, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
void TA2ALowModeCoarseBinned<FImpl, FImplPack, nBasis, binSize>::setup(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &action = envGet(FMat, par().action);
    int         Ls      = env().getObjectLs(par().action);

    if (env().getObjectLs(par().eigenPack) != Ls)
    {
        HADRONS_ERROR(Size, "eigenPack and action Ls mismatch");
    }

    const std::string &schurConv = par().schurConvention;

    if (schurConv == "DiagOne")
    {
        a2a_.reset(new A2AVectorsSchurDiagOne<FImpl>(action));
    }
    else if (schurConv == "DiagTwo")
    {
        a2a_.reset(new A2AVectorsSchurDiagTwo<FImpl>(action));
    }
    else if (schurConv.empty())
    {
        a2a_.reset(new HADRONS_DEFAULT_SCHUR_A2A<FImpl>(action));
    }
    else
    {
        HADRONS_ERROR(Argument, "unknown schurConvention '" + schurConv
                      + "' (expected 'DiagOne', 'DiagTwo', or empty for the compiled-in default)");
    }
    LOG(Message) << "A2A vector construction using Schur convention '"
                 << (schurConv.empty() ? "compiled-in default" : schurConv)
                 << "'" << std::endl;

    if (par().checkInterval > 0)
    {
        LOG(Message) << "In-program eigenvector check enabled: every "
                     << par().checkInterval << " low mode(s), validated "
                     << "against the same operator used to build V/W above, "
                     << "plus a Schur-convention-blind reconstruction check "
                     << "(Mee*v_ie + Meo*v_io == 0, and its W/dagger analogue)"
                     << std::endl;
    }

    Nl_ = epack.evecCoarse.size();
    if (Nl_ % binSize != 0)
    {
        HADRONS_ERROR(Size, "number of low modes (" + std::to_string(Nl_)
                            + ") is not a multiple of the bin size ("
                            + std::to_string(binSize) + ")");
    }
    Nb_ = Nl_ / binSize;

    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(FermionFieldPack, "evecF", Ls, epack.evec[0].Grid());
    envTmp(FermionField, "evecD", Ls, action.FermionRedBlackGrid());
    envTmp(FermionField, "vTmp", 1, envGetGrid(FermionField));
    envTmp(FermionField, "wTmp", 1, envGetGrid(FermionField));
    envTmp(Lattice<SiteSpinorSet>, "vBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
    envTmp(Lattice<SiteSpinorSet>, "wBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
    // Scratch fields for the reconstruction check (see execute()): unlike the
    // eigenvector check above, Mee(v_ie) + Meo(v_io) == 0 is an exact
    // algebraic identity that holds by construction regardless of which
    // Schur convention was used or how good evec_i is as an eigenvector, so
    // it independently validates that makeLowModeV/W's primitive composition
    // is implemented correctly, rather than validating the eigenpack data.
    envTmp(FermionField, "reconCheckE", 1, action.FermionRedBlackGrid());
    envTmp(FermionField, "reconCheckO", 1, action.FermionRedBlackGrid());
    envTmp(FermionField, "reconCheckTmp1", 1, action.FermionRedBlackGrid());
    envTmp(FermionField, "reconCheckTmp2", 1, action.FermionRedBlackGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
void TA2ALowModeCoarseBinned<FImpl, FImplPack, nBasis, binSize>::execute(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &action = envGet(FMat, par().action);
    int         Ls      = env().getObjectLs(par().action);

    envGetTmp(FermionFieldPack, evecF);
    envGetTmp(FermionField, evecD);
    envGetTmp(FermionField, vTmp);
    envGetTmp(FermionField, wTmp);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(Lattice<SiteSpinorSet>, wBin);
    envGetTmp(FermionField, reconCheckE);
    envGetTmp(FermionField, reconCheckO);
    envGetTmp(FermionField, reconCheckTmp1);
    envGetTmp(FermionField, reconCheckTmp2);

    // Mee(v_ie) + Meo(v_io) == 0 (V) and MeeDag(w_ie) + MoeDag(w_io) == 0 (W)
    // are exact algebraic identities that hold by construction for either
    // Schur convention -- see the derivation in memory/conversation history.
    // A residual here isolates a bug in makeLowModeV/W's own primitive
    // composition, independent of the eigenvector check above (which never
    // calls makeLowModeV/W or looks at their output at all) and independent
    // of which Schur convention is in effect (the identity holds either way).
    const RealD reconResidualTol = 1e-10;
    auto checkReconstruction = [&](const FermionField &full, bool dagger,
                                   const std::string &label, unsigned int il)
    {
        pickCheckerboard(Even, reconCheckE, full);
        pickCheckerboard(Odd, reconCheckO, full);
        if (dagger)
        {
            action.MooeeDag(reconCheckE, reconCheckTmp1);
            action.MeooeDag(reconCheckO, reconCheckTmp2);
        }
        else
        {
            action.Mooee(reconCheckE, reconCheckTmp1);
            action.Meooe(reconCheckO, reconCheckTmp2);
        }
        reconCheckTmp1 = reconCheckTmp1 + reconCheckTmp2;

        RealD reconResidual = std::sqrt(norm2(reconCheckTmp1));
        RealD reconScale    = std::sqrt(norm2(reconCheckE));
        RealD reconRelRes   = (reconScale != 0.) ? reconResidual / reconScale : reconResidual;

        LOG(Message) << label << " reconstruction check, low mode " << il
                     << ": |even-block residual| / |even piece| = " << reconRelRes
                     << ((reconRelRes < reconResidualTol) ? "  [OK]" : "  [FAIL]") << std::endl;
    };

    // blockPromote() never sets its output's checkerboard, so evecF would
    // otherwise keep its default (Even) even though it holds the Odd-site
    // data promoted from epack.evec (Odd, see EigenPack.hpp). precisionChange
    // then propagates that wrong tag onto evecD, which corrupts the
    // in-program eigenvector check below (a Schur operator relies on the
    // field's own Checkerboard() for its hopping-term bookkeeping). Does not
    // affect makeLowModeV/W: those force src_o_.Checkerboard() = Odd
    // themselves before use.
    evecF.Checkerboard() = epack.evec[0].Checkerboard();

    LOG(Message) << "Streaming low-mode all-to-all vectors, decompressing/"
                 << "promoting on the fly from compressed coarse eigenpack '"
                 << par().eigenPack << "' (" << Nl_ << " low modes, " << Nb_
                 << " bins of " << binSize << ")" << std::endl;

    for (unsigned int b = 0; b < Nb_; b++)
    {
        for (unsigned int j = 0; j < binSize; j++)
        {
            unsigned int il = b * binSize + j;

            startTimer("Promote");
            blockPromote(epack.evecCoarse[il], evecF, epack.evec);
            precisionChange(evecD, evecF);
            stopTimer("Promote");

            const bool doCheck = (par().checkInterval > 0) && (il % par().checkInterval == 0);

            if (doCheck)
            {
                startTimer("Eigenvector check");

                const RealD                                          checkResidual = 1e-5;
                PlainHermOp<FermionField>                            checkHermOp(a2a_->op());
                ImplicitlyRestartedLanczosHermOpTester<FermionField> checkTester(checkHermOp);
                RealD                                                evalStored = epack.evalCoarse[il];
                RealD                                                evalRecon  = evalStored;
                int conv = checkTester.TestConvergence(il, checkResidual, evecD, evalRecon, 1.0);
                RealD relDiff = (evalStored != 0.)
                              ? std::abs(evalRecon - evalStored) / std::abs(evalStored)
                              : std::abs(evalRecon - evalStored);

                LOG(Message) << "In-program check, low mode " << il
                             << ": stored eval = " << evalStored
                             << ", reconstructed eval = " << evalRecon
                             << ", rel eval diff = " << relDiff
                             << (conv ? "  [OK]" : "  [FAIL]") << std::endl;

                stopTimer("Eigenvector check");
            }

            startTimer("V low mode");
            LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
            if (Ls == 1)
            {
                a2a_->makeLowModeV(vTmp, evecD, epack.evalCoarse[il]);
                if (doCheck)
                {
                    startTimer("V reconstruction check");
                    checkReconstruction(vTmp, false, "V", il);
                    stopTimer("V reconstruction check");
                }
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a_->makeLowModeV5D(vTmp, f5, evecD, epack.evalCoarse[il]);
                if (doCheck)
                {
                    startTimer("V reconstruction check");
                    checkReconstruction(f5, false, "V", il);
                    stopTimer("V reconstruction check");
                }
            }
            pokeLorentz(vBin, vTmp, j);
            stopTimer("V low mode");

            startTimer("W low mode");
            LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
            if (Ls == 1)
            {
                a2a_->makeLowModeW(wTmp, evecD, epack.evalCoarse[il]);
                if (doCheck)
                {
                    startTimer("W reconstruction check");
                    checkReconstruction(wTmp, true, "W", il);
                    stopTimer("W reconstruction check");
                }
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a_->makeLowModeW5D(wTmp, f5, evecD, epack.evalCoarse[il]);
                if (doCheck)
                {
                    startTimer("W reconstruction check");
                    // makeLowModeW5D leaves f5 holding DminusDag(raw w_i), not
                    // the raw w_i itself (unlike makeLowModeV5D, which writes
                    // the raw v_i directly into its vout_5d argument) -- the
                    // reconstruction identity only holds for the raw output,
                    // so re-derive it here with a direct makeLowModeW call.
                    // f5's DminusDag content is no longer needed after
                    // makeLowModeW5D returned (wTmp already holds the real
                    // production output), so it's safe to overwrite.
                    a2a_->makeLowModeW(f5, evecD, epack.evalCoarse[il]);
                    checkReconstruction(f5, true, "W", il);
                    stopTimer("W reconstruction check");
                }
            }
            pokeLorentz(wBin, wTmp, j);
            stopTimer("W low mode");
        }

        startTimer("V I/O");
        A2AVectorsIo::writeElement(par().output + "_v", vBin, b, vm().getTrajectory());
        stopTimer("V I/O");

        startTimer("W I/O");
        A2AVectorsIo::writeElement(par().output + "_w", wBin, b, vm().getTrajectory());
        stopTimer("W I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ALowModeCoarseBinned_hpp_
