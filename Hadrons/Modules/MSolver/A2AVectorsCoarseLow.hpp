/*
 * A2AVectorsCoarseLow.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MSolver_A2AVectorsCoarseLow_hpp_
#define Hadrons_MSolver_A2AVectorsCoarseLow_hpp_

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
 *  from a compressed (LC) coarse eigenpack. Low modes are decompressed/     *
 *  promoted on the fly, one bin at a time, so the full Nl-sized V/W arrays  *
 *  are never resident -- unlike MSolver::A2AVectorsCoarse followed by       *
 *  MIO::SaveBinnedA2AVecs, which together hold the full arrays twice (once  *
 *  as the flat V/W vectors, again as the packed bin copy). Only one bin's   *
 *  worth of fields (compile-time binSize) is ever resident at once.         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsCoarseLowPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsCoarseLowPar,
                                    std::string, eigenPack,
                                    std::string, action,
                                    std::string, solver,
                                    std::string, output,
                                    std::string, schurConvention,
                                    unsigned int, checkInterval);
};

template <typename FImpl, typename FImplPack, int nBasis, int binSize>
class TA2AVectorsCoarseLow : public Module<A2AVectorsCoarseLowPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    FERM_TYPE_ALIASES(FImplPack, Pack);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl>          A2A;
    typedef CoarseFermionEigenPack<FImplPack, nBasis> EPack;
    typedef typename FImpl::SiteSpinor::vector_type   vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2AVectorsCoarseLow(const std::string name);
    // destructor
    virtual ~TA2AVectorsCoarseLow(void) {};
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
    // built once in setup() from 'action' + par().schurConvention when
    // par().checkInterval > 0; held here (rather than rebuilt per mode) since
    // neither the action nor the convention change over the module's
    // lifetime. See MSolver::MixedPrecisionRBPrecCG for why this has to be a
    // pointer to the shared base class rather than a concrete stack object:
    // the concrete Schur convention is only known at runtime (from XML), and
    // a C++ variable's type must be fixed at compile time.
    std::unique_ptr<SchurOperatorBase<FermionField>> checkOp_;
};

MODULE_REGISTER_TMP(A2AVectorsCoarseLow200Bin200,
    ARG(TA2AVectorsCoarseLow<FIMPL, FIMPLF, 200, 200>), MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseLow200Bin100,
    ARG(TA2AVectorsCoarseLow<FIMPL, FIMPLF, 200, 100>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseLow200Bin200,
    ARG(TA2AVectorsCoarseLow<ZFIMPL, ZFIMPLF, 200, 200>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseLow200Bin100,
    ARG(TA2AVectorsCoarseLow<ZFIMPL, ZFIMPLF, 200, 100>), MSolver);

/******************************************************************************
 *                   TA2AVectorsCoarseLow implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
TA2AVectorsCoarseLow<FImpl, FImplPack, nBasis, binSize>::TA2AVectorsCoarseLow(const std::string name)
: Module<A2AVectorsCoarseLowPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
std::vector<std::string> TA2AVectorsCoarseLow<FImpl, FImplPack, nBasis, binSize>::getInput(void)
{
    // Unlike MSolver::A2AVectorsCoarse, this module never runs a CG solve
    // (low modes are exact deflation, not stochastic), so par().solver is
    // only referenced to satisfy A2AVectorsSchurDiagTwo's constructor -- it
    // is never actually invoked. Point it at whatever solver object is
    // otherwise convenient; there is no "_subtract" requirement here.
    std::vector<std::string> in = {par().eigenPack, par().action, par().solver};

    return in;
}

template <typename FImpl, typename FImplPack, int nBasis, int binSize>
std::vector<std::string> TA2AVectorsCoarseLow<FImpl, FImplPack, nBasis, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
void TA2AVectorsCoarseLow<FImpl, FImplPack, nBasis, binSize>::setup(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &action = envGet(FMat, par().action);
    auto        &solver = envGet(Solver, par().solver);
    int         Ls      = env().getObjectLs(par().action);

    if (env().getObjectLs(par().eigenPack) != Ls)
    {
        HADRONS_ERROR(Size, "eigenPack and action Ls mismatch");
    }

    if (par().checkInterval > 0)
    {
        const std::string &schurConv = par().schurConvention;

        if (schurConv == "DiagOne")
        {
            checkOp_.reset(new SchurDiagOneOperator<FMat, FermionField>(action));
        }
        else if (schurConv == "DiagTwo")
        {
            checkOp_.reset(new SchurDiagTwoOperator<FMat, FermionField>(action));
        }
        else if (schurConv.empty())
        {
            checkOp_.reset(new HADRONS_DEFAULT_SCHUR_OP<FMat, FermionField>(action));
        }
        else
        {
            HADRONS_ERROR(Argument, "unknown schurConvention '" + schurConv
                          + "' (expected 'DiagOne', 'DiagTwo', or empty for the compiled-in default)");
        }
        LOG(Message) << "In-program eigenvector check enabled: every "
                     << par().checkInterval << " low mode(s), Schur convention '"
                     << (schurConv.empty() ? "compiled-in default" : schurConv)
                     << "'" << std::endl;
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
    envTmp(A2A, "a2a", 1, action, solver);
    envTmp(FermionFieldPack, "evecF", Ls, epack.evec[0].Grid());
    envTmp(FermionField, "evecD", Ls, action.FermionRedBlackGrid());
    envTmp(FermionField, "vTmp", 1, envGetGrid(FermionField));
    envTmp(FermionField, "wTmp", 1, envGetGrid(FermionField));
    envTmp(Lattice<SiteSpinorSet>, "vBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
    envTmp(Lattice<SiteSpinorSet>, "wBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis, int binSize>
void TA2AVectorsCoarseLow<FImpl, FImplPack, nBasis, binSize>::execute(void)
{
    auto        &epack = envGet(EPack, par().eigenPack);
    int         Ls     = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);
    envGetTmp(FermionFieldPack, evecF);
    envGetTmp(FermionField, evecD);
    envGetTmp(FermionField, vTmp);
    envGetTmp(FermionField, wTmp);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(Lattice<SiteSpinorSet>, wBin);

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

            if ((par().checkInterval > 0) && (il % par().checkInterval == 0))
            {
                startTimer("Eigenvector check");

                const RealD                                          checkResidual = 1e-5;
                PlainHermOp<FermionField>                            checkHermOp(*checkOp_);
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
                a2a.makeLowModeV(vTmp, evecD, epack.evalCoarse[il]);
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a.makeLowModeV5D(vTmp, f5, evecD, epack.evalCoarse[il]);
            }
            pokeLorentz(vBin, vTmp, j);
            stopTimer("V low mode");

            startTimer("W low mode");
            LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
            if (Ls == 1)
            {
                a2a.makeLowModeW(wTmp, evecD, epack.evalCoarse[il]);
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a.makeLowModeW5D(wTmp, f5, evecD, epack.evalCoarse[il]);
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

#endif // Hadrons_MSolver_A2AVectorsCoarseLow_hpp_
