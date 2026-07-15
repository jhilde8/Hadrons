/*
 * PrecisionCastCoarseEPack.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
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
#ifndef Hadrons_MUtilities_PrecisionCastCoarseEPack_hpp_
#define Hadrons_MUtilities_PrecisionCastCoarseEPack_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *              Precision cast module for compressed (LC) eigenpacks         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class PrecisionCastCoarseEPackPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PrecisionCastCoarseEPackPar,
                                    std::string, in,
                                    std::string, blockSize,
                                    bool,        redBlack);
};

template <typename FImplIn, typename FImplOut, int nBasis>
class TPrecisionCastCoarseEPack: public Module<PrecisionCastCoarseEPackPar>
{
public:
    typedef CoarseFermionEigenPack<FImplIn,  nBasis> PackIn;
    typedef CoarseFermionEigenPack<FImplOut, nBasis> PackOut;
    typedef typename PackOut::Field                  FieldOut;
    typedef typename PackOut::FieldIo                FieldIoOut;
    typedef typename PackOut::CoarseField            CoarseFieldOut;
    typedef typename PackOut::CoarseFieldIo          CoarseFieldIoOut;
    typedef CoarseEigenPack<FieldOut, CoarseFieldOut, FieldIoOut, CoarseFieldIoOut> BasePack;

public:
    // constructor
    TPrecisionCastCoarseEPack(const std::string name);
    // destructor
    virtual ~TPrecisionCastCoarseEPack(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CoarseFermionDoublePrecisionCastEPack200,
                    ARG(TPrecisionCastCoarseEPack<FIMPLF, FIMPLD, 200>),
                    MUtilities);

/******************************************************************************
 *                 TPrecisionCastCoarseEPack implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut, int nBasis>
TPrecisionCastCoarseEPack<FImplIn, FImplOut, nBasis>::TPrecisionCastCoarseEPack(const std::string name)
: Module<PrecisionCastCoarseEPackPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut, int nBasis>
std::vector<std::string> TPrecisionCastCoarseEPack<FImplIn, FImplOut, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().in};

    return in;
}

template <typename FImplIn, typename FImplOut, int nBasis>
std::vector<std::string> TPrecisionCastCoarseEPack<FImplIn, FImplOut, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut, int nBasis>
void TPrecisionCastCoarseEPack<FImplIn, FImplOut, nBasis>::setup(void)
{
    auto         &in        = envGet(PackIn, par().in);
    auto         blockSize  = strToVec<int>(par().blockSize);
    unsigned int Ls         = env().getObjectLs(par().in);
    GridBase     *gridIo       = nullptr;
    GridBase     *gridCoarseIo = nullptr;
    GridBase     *fineGrid   = getGrid<FieldOut>(par().redBlack, Ls);
    GridBase     *coarseGrid = envGetCoarseGrid(CoarseFieldOut, blockSize, Ls);

    envCreateDerived(BasePack, PackOut, getName(), Ls,
                     in.evec.size(), in.evecCoarse.size(),
                     fineGrid, coarseGrid, gridIo, gridCoarseIo);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImplIn, typename FImplOut, int nBasis>
void TPrecisionCastCoarseEPack<FImplIn, FImplOut, nBasis>::execute(void)
{
    LOG(Message) << "Casting compressed eigenpack '" << par().in << "'" << std::endl;
    LOG(Message) << "In  type: " << typeName<PackIn>()  << std::endl;
    LOG(Message) << "Out type: " << typeName<PackOut>() << std::endl;

    auto &in  = envGet(PackIn,  par().in);
    auto &out = envGetDerived(BasePack, PackOut, getName());

    out.record = in.record;

    LOG(Message) << "Casting " << in.evec.size() << " fine basis vector(s)" << std::endl;
    for (unsigned int i = 0; i < in.evec.size(); ++i)
    {
        precisionChange(out.evec[i], in.evec[i]);
        out.eval[i] = in.eval[i];
    }

    LOG(Message) << "Casting " << in.evecCoarse.size() << " coarse coefficient vector(s)" << std::endl;
    for (unsigned int i = 0; i < in.evecCoarse.size(); ++i)
    {
        precisionChange(out.evecCoarse[i], in.evecCoarse[i]);
        out.evalCoarse[i] = in.evalCoarse[i];
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_PrecisionCastCoarseEPack_hpp_
