/*
 * A2AVectorsCoarse.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MSolver_A2AVectorsCoarse_hpp_
#define Hadrons_MSolver_A2AVectorsCoarse_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Create all-to-all V & W vectors directly from a compressed (LC) coarse   *
 *  eigenpack -- low modes are decompressed/promoted on the fly, one at a    *
 *  time, instead of reading a resident dense eigenvector array. High modes  *
 *  are unchanged from MSolver::A2AVectors (noise and the mixed-precision    *
 *  solver output are already double precision).                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsCoarsePar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsCoarsePar,
                                  std::string, noise,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, solver,
                                  std::string, output,
                                  bool,        multiFile);
};

template <typename FImpl, typename FImplPack, int nBasis>
class TA2AVectorsCoarse : public Module<A2AVectorsCoarsePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    FERM_TYPE_ALIASES(FImplPack, Pack);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl>          A2A;
    typedef CoarseFermionEigenPack<FImplPack, nBasis> EPack;
public:
    // constructor
    TA2AVectorsCoarse(const std::string name);
    // destructor
    virtual ~TA2AVectorsCoarse(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
};

MODULE_REGISTER_TMP(A2AVectorsCoarse200,
    ARG(TA2AVectorsCoarse<FIMPL, FIMPLF, 200>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarse200,
    ARG(TA2AVectorsCoarse<ZFIMPL, ZFIMPLF, 200>), MSolver);

/******************************************************************************
 *                     TA2AVectorsCoarse implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
TA2AVectorsCoarse<FImpl, FImplPack, nBasis>::TA2AVectorsCoarse(const std::string name)
: Module<A2AVectorsCoarsePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
std::vector<std::string> TA2AVectorsCoarse<FImpl, FImplPack, nBasis>::getInput(void)
{
    // Always has low modes (that is the point of this module), so -- same as
    // MSolver::A2AVectors's hasLowModes branch -- the solver used is always
    // the "_subtract" variant (deflation contribution already subtracted
    // out, so the high-mode CG isn't redundant with the low-mode
    // construction below).
    std::vector<std::string> in = {par().eigenPack, par().solver + "_subtract",
                                   par().noise};

    return in;
}

template <typename FImpl, typename FImplPack, int nBasis>
std::vector<std::string> TA2AVectorsCoarse<FImpl, FImplPack, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName() + "_v", getName() + "_w"};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
void TA2AVectorsCoarse<FImpl, FImplPack, nBasis>::setup(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &noise  = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action = envGet(FMat, par().action);
    auto        &solver = envGet(Solver, par().solver + "_subtract");
    int         Ls      = env().getObjectLs(par().action);

    if (env().getObjectLs(par().eigenPack) != Ls)
    {
        HADRONS_ERROR(Size, "eigenPack and action Ls mismatch");
    }

    Nl_ = epack.evecCoarse.size();

    envCreate(std::vector<FermionField>, getName() + "_v", 1,
              Nl_ + noise.fermSize(), envGetGrid(FermionField));
    envCreate(std::vector<FermionField>, getName() + "_w", 1,
              Nl_ + noise.fermSize(), envGetGrid(FermionField));

    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
    envTmp(FermionFieldPack, "evecF", Ls, epack.evec[0].Grid());
    envTmp(FermionField, "evecD", Ls, action.FermionRedBlackGrid());
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
void TA2AVectorsCoarse<FImpl, FImplPack, nBasis>::execute(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &noise  = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &v      = envGet(std::vector<FermionField>, getName() + "_v");
    auto        &w      = envGet(std::vector<FermionField>, getName() + "_w");
    int         Ls      = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);
    envGetTmp(FermionFieldPack, evecF);
    envGetTmp(FermionField, evecD);

    LOG(Message) << "Computing all-to-all vectors, decompressing/promoting "
                 << "low modes on the fly from compressed coarse eigenpack '"
                 << par().eigenPack << "' (" << Nl_ << " low modes) and "
                 << "noise '" << par().noise << "' (" << noise.fermSize()
                 << " noise vectors)" << std::endl;

    // Low modes: promote each stored eigenvector from the compressed coarse
    // representation on the fly -- single precision, forced by the shared
    // subspace -- cast the one promoted vector up to double, then feed the
    // existing double-precision low-mode construction unchanged.
    for (unsigned int il = 0; il < Nl_; il++)
    {
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        blockPromote(epack.evecCoarse[il], evecF, epack.evec);
        precisionChange(evecD, evecF);
        if (Ls == 1)
        {
            a2a.makeLowModeV(v[il], evecD, epack.evalCoarse[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v[il], f5, evecD, epack.evalCoarse[il]);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w[il], evecD, epack.evalCoarse[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w[il], f5, evecD, epack.evalCoarse[il]);
        }
        stopTimer("W low mode");
    }

    // High modes: unchanged from MSolver::A2AVectors.
    for (unsigned int ih = 0; ih < noise.fermSize(); ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "")
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("V high mode");
        startTimer("W high mode");
        LOG(Message) << "W vector i = " << Nl_ + ih
                     << " (" << ((Nl_ > 0) ? "high " : "")
                     << "stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeW(w[Nl_ + ih], noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeW5D(w[Nl_ + ih], f5, noise.getFerm(ih));
        }
        stopTimer("W high mode");
    }

    // I/O if necessary
    if (!par().output.empty())
    {
        startTimer("V I/O");
        A2AVectorsIo::write(par().output + "_v", v, par().multiFile, vm().getTrajectory());
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::write(par().output + "_w", w, par().multiFile, vm().getTrajectory());
        stopTimer("W I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsCoarse_hpp_
