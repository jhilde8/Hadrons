/*
 * A2AVectorsExactLow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MSolver_A2AVectorsExactLow_hpp_
#define Hadrons_MSolver_A2AVectorsExactLow_hpp_

#include <memory>

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Stream low-mode A2A V & W vectors straight to disk, binned, directly     *
 *  from a resident, single-precision exact (uncompressed) eigenpack --      *
 *  companion to MSolver::A2AVectorsCoarseLow for the exact-deflation case.  *
 *  Unlike a coarse/LC-compressed eigenpack there is no basis to decompress  *
 *  from: each low mode only needs a plain precisionChange() promotion to    *
 *  double, one mode at a time, so the full Nl-sized double-precision array  *
 *  is never resident -- only the single-precision eigenpack (already kept   *
 *  around for MGuesser::ExactDeflationF's inner-CG guesser) stays in        *
 *  memory, plus one bin's worth (compile-time binSize) of promoted fields.  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsExactLowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsExactLowPar,
                                    std::string, eigenPack,
                                    std::string, action,
                                    std::string, solver,
                                    std::string, output,
                                    std::string, schurConvention,
                                    unsigned int, checkInterval);
};

template <typename FImpl, typename FImplPack, int binSize>
class TA2AVectorsExactLow : public Module<A2AVectorsExactLowPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef BaseFermionEigenPack<FImplPack>         EPack;
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2AVectorsExactLow(const std::string name);
    // destructor
    virtual ~TA2AVectorsExactLow(void) {};
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

MODULE_REGISTER_TMP(A2AVectorsExactLowBin200,
    ARG(TA2AVectorsExactLow<FIMPL, FIMPLF, 200>), MSolver);
MODULE_REGISTER_TMP(A2AVectorsExactLowBin100,
    ARG(TA2AVectorsExactLow<FIMPL, FIMPLF, 100>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsExactLowBin200,
    ARG(TA2AVectorsExactLow<ZFIMPL, ZFIMPLF, 200>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsExactLowBin100,
    ARG(TA2AVectorsExactLow<ZFIMPL, ZFIMPLF, 100>), MSolver);

/******************************************************************************
 *                   TA2AVectorsExactLow implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int binSize>
TA2AVectorsExactLow<FImpl, FImplPack, binSize>::TA2AVectorsExactLow(const std::string name)
: Module<A2AVectorsExactLowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int binSize>
std::vector<std::string> TA2AVectorsExactLow<FImpl, FImplPack, binSize>::getInput(void)
{
    // Unlike MSolver::A2AVectorsCoarse, this module never runs a CG solve
    // (low modes are exact deflation, not stochastic), so par().solver is
    // only referenced to satisfy the A2A vector constructor -- it is never
    // actually invoked. Point it at whatever solver object is otherwise
    // convenient; there is no "_subtract" requirement here.
    std::vector<std::string> in = {par().eigenPack, par().action, par().solver};

    return in;
}

template <typename FImpl, typename FImplPack, int binSize>
std::vector<std::string> TA2AVectorsExactLow<FImpl, FImplPack, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int binSize>
void TA2AVectorsExactLow<FImpl, FImplPack, binSize>::setup(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &action = envGet(FMat, par().action);
    auto        &solver = envGet(Solver, par().solver);
    int         Ls      = env().getObjectLs(par().action);

    if (env().getObjectLs(par().eigenPack) != Ls)
    {
        HADRONS_ERROR(Size, "eigenPack and action Ls mismatch");
    }

    const std::string &schurConv = par().schurConvention;

    if (schurConv == "DiagOne")
    {
        a2a_.reset(new A2AVectorsSchurDiagOne<FImpl>(action, solver));
    }
    else if (schurConv == "DiagTwo")
    {
        a2a_.reset(new A2AVectorsSchurDiagTwo<FImpl>(action, solver));
    }
    else if (schurConv.empty())
    {
        a2a_.reset(new HADRONS_DEFAULT_SCHUR_A2A<FImpl>(action, solver));
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
                     << "against the same operator used to build V/W above"
                     << std::endl;
    }

    Nl_ = epack.evec.size();
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
    envTmp(FermionField, "evecD", Ls, action.FermionRedBlackGrid());
    envTmp(FermionField, "vTmp", 1, envGetGrid(FermionField));
    envTmp(FermionField, "wTmp", 1, envGetGrid(FermionField));
    envTmp(Lattice<SiteSpinorSet>, "vBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
    envTmp(Lattice<SiteSpinorSet>, "wBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int binSize>
void TA2AVectorsExactLow<FImpl, FImplPack, binSize>::execute(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &action = envGet(FMat, par().action);
    int         Ls      = env().getObjectLs(par().action);

    envGetTmp(FermionField, evecD);
    envGetTmp(FermionField, vTmp);
    envGetTmp(FermionField, wTmp);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(Lattice<SiteSpinorSet>, wBin);

    LOG(Message) << "Streaming low-mode all-to-all vectors, promoting on the "
                 << "fly from resident single-precision exact eigenpack '"
                 << par().eigenPack << "' (" << Nl_ << " low modes, " << Nb_
                 << " bins of " << binSize << ")" << std::endl;

    for (unsigned int b = 0; b < Nb_; b++)
    {
        for (unsigned int j = 0; j < binSize; j++)
        {
            unsigned int il = b * binSize + j;

            startTimer("Promote");
            precisionChange(evecD, epack.evec[il]);
            stopTimer("Promote");

            const bool doCheck = (par().checkInterval > 0) && (il % par().checkInterval == 0);

            if (doCheck)
            {
                startTimer("Eigenvector check");

                const RealD                                          checkResidual = 1e-5;
                PlainHermOp<FermionField>                            checkHermOp(a2a_->op());
                ImplicitlyRestartedLanczosHermOpTester<FermionField> checkTester(checkHermOp);
                RealD                                                evalStored = epack.eval[il];
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
                a2a_->makeLowModeV(vTmp, evecD, epack.eval[il]);
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a_->makeLowModeV5D(vTmp, f5, evecD, epack.eval[il]);
            }
            pokeLorentz(vBin, vTmp, j);
            stopTimer("V low mode");

            startTimer("W low mode");
            LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
            if (Ls == 1)
            {
                a2a_->makeLowModeW(wTmp, evecD, epack.eval[il]);
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a_->makeLowModeW5D(wTmp, f5, evecD, epack.eval[il]);
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

#endif // Hadrons_MSolver_A2AVectorsExactLow_hpp_
