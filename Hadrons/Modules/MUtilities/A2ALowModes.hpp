/*
 * A2ALowModes.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef Hadrons_MUtilities_A2ALowModes_hpp_
#define Hadrons_MUtilities_A2ALowModes_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                  Create low-mode all-to-all V & W vectors                  *
 ******************************************************************************/
/*
 * Each low mode's V and W are written as soon as they are built, one fermion
 * field per file at <output>_v.<traj>/elem<i>.bin and <output>_w.<traj>/
 * elem<i>.bin, i the eigenpack index. An empty output computes without
 * writing.
 */
BEGIN_MODULE_NAMESPACE(MUtilities)
class A2ALowModesPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALowModesPar,
                                  std::string, action,
                                  std::string, eigenPack,
                                  std::string, output);
};
template <typename FImpl, typename Pack>
class TA2ALowModes : public Module<A2ALowModesPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2ALowModes(const std::string name);
    // destructor
    virtual ~TA2ALowModes(void) {};
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
MODULE_REGISTER_TMP(A2ALowModes,
    ARG(TA2ALowModes<FIMPL, BaseFermionEigenPack<FIMPL>>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModes,
    ARG(TA2ALowModes<ZFIMPL, BaseFermionEigenPack<ZFIMPL>>), MUtilities);
/******************************************************************************
 *                       TA2ALowModes implementation                          *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
TA2ALowModes<FImpl, Pack>::TA2ALowModes(const std::string name)
: Module<A2ALowModesPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
std::vector<std::string> TA2ALowModes<FImpl, Pack>::getInput(void)
{
    assert(!par().eigenPack.empty());
    std::vector<std::string> in = {par().action, par().eigenPack};

    return in;
}

template <typename FImpl, typename Pack>
std::vector<std::string> TA2ALowModes<FImpl, Pack>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2ALowModes<FImpl, Pack>::setup(void)
{
    assert(!par().eigenPack.empty());
    auto &action = envGet(FMat, par().action);
    auto &epack  = envGet(Pack, par().eigenPack);
    int  Ls      = env().getObjectLs(par().action);

    Nl_ = epack.evec.size();
    envTmpLat(FermionField, "v");
    envTmpLat(FermionField, "w");
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename Pack>
void TA2ALowModes<FImpl, Pack>::execute(void)
{
    int Ls = env().getObjectLs(par().action);
    envGetTmp(FermionField, v);
    envGetTmp(FermionField, w);
    envGetTmp(A2A, a2a);

    LOG(Message) << "Getting low-mode part of all-to-all vectors "
                 << " from '" << par().eigenPack << "' ("
                 << Nl_ << " eigenPack modes)" << std::endl;

    // Low modes
    for (unsigned int il = 0; il < Nl_; il++)
    {
        auto &epack = envGet(Pack, par().eigenPack);
        startTimer("V low mode");
        LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeV(v, epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeV5D(v, f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("V low mode");
        startTimer("W low mode");
        LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeLowModeW(w, epack.evec[il], epack.eval[il]);
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeLowModeW5D(w, f5, epack.evec[il], epack.eval[il]);
        }
        stopTimer("W low mode");

        if (!par().output.empty())
        {
            startTimer("V I/O");
            A2AVectorsIo::writeElement(par().output + "_v", v, il,
                                       vm().getTrajectory());
            stopTimer("V I/O");
            startTimer("W I/O");
            A2AVectorsIo::writeElement(par().output + "_w", w, il,
                                       vm().getTrajectory());
            stopTimer("W I/O");
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MUtilities_A2ALowModes_hpp_
