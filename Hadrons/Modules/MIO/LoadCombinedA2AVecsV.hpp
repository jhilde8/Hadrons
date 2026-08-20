/*
 * LoadCombinedA2AVecsV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MIO_LoadCombinedA2AVecsV_hpp_
#define Hadrons_MIO_LoadCombinedA2AVecsV_hpp_

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

class LoadCombinedA2AVecsVPar: Serializable
{
public:
    // nHit: 1/nHit hit-average normalization applied to the high-mode blocks
    // after loading, matching MIO::LoadBinnedA2AVecsV -- the multi-hit
    // estimator averages over hits, and by convention the whole factor lives
    // on the V side while W stays raw unit-modulus noise. Set it to the hit
    // count only when loading V vectors; leave it 0 (absent from the XML) or
    // 1 -- both mean no rescaling -- for W loads and single-hit data.
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCombinedA2AVecsVPar,
                                    std::string,               lowFilestem,
                                    unsigned int,              nLow,
                                    std::string,               highStem,
                                    std::vector<std::string>, highExtensions,
                                    unsigned int,              highSize,
                                    unsigned int,              nHit);
};

template <typename FImpl, int lowBinSize, int highBinSize>
class TLoadCombinedA2AVecsV: public Module<LoadCombinedA2AVecsVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef iVector<typename FImpl::SiteSpinor, lowBinSize>  LowBinnedSpinor;
    typedef iVector<typename FImpl::SiteSpinor, highBinSize> HighBinnedSpinor;
public:
    // constructor
    TLoadCombinedA2AVecsV(const std::string name);
    // destructor
    virtual ~TLoadCombinedA2AVecsV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCombinedA2AVecsV200x128, ARG(TLoadCombinedA2AVecsV<FIMPL, 200, 128>), MIO);
MODULE_REGISTER_TMP(LoadCombinedA2AVecsV100x128, ARG(TLoadCombinedA2AVecsV<FIMPL, 100, 128>), MIO);
MODULE_REGISTER_TMP(LoadCombinedA2AVecsV200x96,  ARG(TLoadCombinedA2AVecsV<FIMPL, 200, 96>),  MIO);
MODULE_REGISTER_TMP(LoadCombinedA2AVecsV100x96,  ARG(TLoadCombinedA2AVecsV<FIMPL, 100, 96>),  MIO);

/******************************************************************************
 *                  TLoadCombinedA2AVecsV implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
TLoadCombinedA2AVecsV<FImpl, lowBinSize, highBinSize>::TLoadCombinedA2AVecsV(const std::string name)
: Module<LoadCombinedA2AVecsVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
std::vector<std::string> TLoadCombinedA2AVecsV<FImpl, lowBinSize, highBinSize>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl, int lowBinSize, int highBinSize>
std::vector<std::string> TLoadCombinedA2AVecsV<FImpl, lowBinSize, highBinSize>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
void TLoadCombinedA2AVecsV<FImpl, lowBinSize, highBinSize>::setup(void)
{
    if (par().nLow % lowBinSize != 0)
    {
        HADRONS_ERROR(Size, "nLow (" + std::to_string(par().nLow)
                            + ") is not a multiple of the low-mode bin size ("
                            + std::to_string(lowBinSize) + ")");
    }
    if (par().highSize % highBinSize != 0)
    {
        HADRONS_ERROR(Size, "highSize (" + std::to_string(par().highSize)
                            + ") is not a multiple of the high-mode bin size ("
                            + std::to_string(highBinSize) + ")");
    }

    unsigned int total = par().nLow + par().highExtensions.size() * par().highSize;
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
void TLoadCombinedA2AVecsV<FImpl, lowBinSize, highBinSize>::execute(void)
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

    int  Nb   = par().highSize / highBinSize;
    Real norm = (par().nHit > 0) ? (1.0/par().nHit) : 1.0;

    for (unsigned int h = 0; h < par().highExtensions.size(); ++h)
    {
        std::string filestem = par().highStem + par().highExtensions[h];
        std::vector<Lattice<HighBinnedSpinor>> bvec(Nb, envGetGrid(Lattice<HighBinnedSpinor>));

        LOG(Message) << "Loading " << par().highSize << " high-mode A2A vectors from '"
                     << filestem << "' (" << Nb << " files of "
                     << highBinSize << " binned vectors)" << std::endl;
        A2AVectorsIo::read(bvec, filestem, true, vm().getTrajectory());
        A2Autils<FImpl>::template UnpackBinnedVectors<highBinSize>(out, offset, bvec);
        if (norm != 1.0)
        {
            LOG(Message) << "Applying 1/nHit = " << norm
                         << " hit-average normalization" << std::endl;
            for (unsigned int i = offset; i < offset + par().highSize; ++i)
            {
                out[i] = norm*out[i];
            }
        }
        offset += par().highSize;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCombinedA2AVecsV_hpp_
