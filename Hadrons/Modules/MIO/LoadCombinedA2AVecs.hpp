/*
 * LoadCombinedA2AVecs.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
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
#ifndef Hadrons_MIO_LoadCombinedA2AVecs_hpp_
#define Hadrons_MIO_LoadCombinedA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Load and directly assemble a full A2A vector array: an optional block of *
 *  low modes (nLow > 0) followed by one or more blocks of high modes, one   *
 *  per hit -- e.g. light gets low modes + N hits of high modes; strange/    *
 *  charm (no low modes, nLow = 0, lowFilestem unused) get just N hits of    *
 *  high modes concatenated. Every vector is read from disk straight into    *
 *  its final resting slot of one array allocated up front -- unlike         *
 *  Load+Load+...+MUtilities::CombineA2AVecs, there are no per-source        *
 *  temporary arrays and no extra copy pass to splice them together after.  *
 *                                                                            *
 *  High-mode files are always assumed multi-file (per-vector), since a      *
 *  single-file high-mode set on this ensemble's volume is not a realistic   *
 *  case this module needs to worry about. The per-hit high filestem is      *
 *  'highStem + highExtensions[h]' -- e.g. highStem = ".../vw/" and          *
 *  highExtensions = {"l0_v", "l1_v", ...} -- so the common directory only   *
 *  has to be written once rather than once per hit.                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCombinedA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCombinedA2AVecsPar,
                                    std::string,               lowFilestem,
                                    unsigned int,              nLow,
                                    std::string,               highStem,
                                    std::vector<std::string>, highExtensions,
                                    unsigned int,              nHighEach);
};

template <typename FImpl, int lowBinSize, int highBinSize>
class TLoadCombinedA2AVecs: public Module<LoadCombinedA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef iVector<typename FImpl::SiteSpinor, lowBinSize>  LowBinnedSpinor;
    typedef iVector<typename FImpl::SiteSpinor, highBinSize> HighBinnedSpinor;
public:
    // constructor
    TLoadCombinedA2AVecs(const std::string name);
    // destructor
    virtual ~TLoadCombinedA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCombinedA2AVecs200x128, ARG(TLoadCombinedA2AVecs<FIMPL, 200, 128>), MIO);
MODULE_REGISTER_TMP(LoadCombinedA2AVecs100x128, ARG(TLoadCombinedA2AVecs<FIMPL, 100, 128>), MIO);
MODULE_REGISTER_TMP(LoadCombinedA2AVecs200x96,  ARG(TLoadCombinedA2AVecs<FIMPL, 200, 96>),  MIO);
MODULE_REGISTER_TMP(LoadCombinedA2AVecs100x96,  ARG(TLoadCombinedA2AVecs<FIMPL, 100, 96>),  MIO);

/******************************************************************************
 *                  TLoadCombinedA2AVecs implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
TLoadCombinedA2AVecs<FImpl, lowBinSize, highBinSize>::TLoadCombinedA2AVecs(const std::string name)
: Module<LoadCombinedA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
std::vector<std::string> TLoadCombinedA2AVecs<FImpl, lowBinSize, highBinSize>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl, int lowBinSize, int highBinSize>
std::vector<std::string> TLoadCombinedA2AVecs<FImpl, lowBinSize, highBinSize>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
void TLoadCombinedA2AVecs<FImpl, lowBinSize, highBinSize>::setup(void)
{
    if (par().nLow % lowBinSize != 0)
    {
        HADRONS_ERROR(Size, "nLow (" + std::to_string(par().nLow)
                            + ") is not a multiple of the low-mode bin size ("
                            + std::to_string(lowBinSize) + ")");
    }
    if (par().nHighEach % highBinSize != 0)
    {
        HADRONS_ERROR(Size, "nHighEach (" + std::to_string(par().nHighEach)
                            + ") is not a multiple of the high-mode bin size ("
                            + std::to_string(highBinSize) + ")");
    }

    unsigned int total = par().nLow + par().highExtensions.size() * par().nHighEach;
    auto         grid  = envGetGrid(FermionField);

    envCreate(std::vector<FermionField>, getName(), 1, 0, grid);
    auto &out = envGet(std::vector<FermionField>, getName());
    out.reserve(total);
    for (unsigned int i = 0; i < total; ++i)
    {
        out.emplace_back(grid, CpuWrite);
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
void TLoadCombinedA2AVecs<FImpl, lowBinSize, highBinSize>::execute(void)
{
    auto &out = envGet(std::vector<FermionField>, getName());
    unsigned int offset = 0;

    if (par().nLow > 0)
    {
        int Nb = par().nLow / lowBinSize;
        std::vector<Lattice<LowBinnedSpinor>> bvec(Nb, envGetGrid(Lattice<LowBinnedSpinor>));

        LOG(Message) << "Loading " << par().nLow << " low-mode A2A vectors from "
                     << Nb << " files of " << lowBinSize << " binned vectors" << std::endl;
        A2AVectorsIo::read(bvec, par().lowFilestem, true, vm().getTrajectory());
        A2Autils<FImpl>::template UnpackBinnedVectors<lowBinSize>(out, offset, bvec);
        offset += par().nLow;
    }

    int Nb = par().nHighEach / highBinSize;
    for (unsigned int h = 0; h < par().highExtensions.size(); ++h)
    {
        std::string filestem = par().highStem + par().highExtensions[h];
        std::vector<Lattice<HighBinnedSpinor>> bvec(Nb, envGetGrid(Lattice<HighBinnedSpinor>));

        LOG(Message) << "Loading " << par().nHighEach << " high-mode A2A vectors from '"
                     << filestem << "' (" << Nb << " files of "
                     << highBinSize << " binned vectors)" << std::endl;
        A2AVectorsIo::read(bvec, filestem, true, vm().getTrajectory());
        A2Autils<FImpl>::template UnpackBinnedVectors<highBinSize>(out, offset, bvec);
        offset += par().nHighEach;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCombinedA2AVecs_hpp_
