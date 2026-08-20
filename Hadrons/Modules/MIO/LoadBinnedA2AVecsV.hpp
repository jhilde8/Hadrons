/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MIO/LoadBinnedA2AVecsV.hpp

Copyright (C) 2015-2019
Author: Antonin Portelli <antonin.portelli@me.com>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MIO_LoadBinnedA2AVecsV_hpp_
#define Hadrons_MIO_LoadBinnedA2AVecsV_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                 Module to load all-to-all V vectors                        *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)
class LoadBinnedA2AVecsVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadBinnedA2AVecsVPar,
                                    std::string,              lowFilestem,
                                    std::vector<std::string>, highFileStems,
                                    bool,                     multiFile,
                                    unsigned int,             lowSize,
                                    unsigned int,             highSize,
                                    unsigned int,             nHit);
};
template <typename FImpl, int lowBinSize, int highBinSize>
class TLoadBinnedA2AVecsV: public Module<LoadBinnedA2AVecsVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    // typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, lowBinSize>
    //     LowSiteSpinorSet;
    // typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, highBinSize>
    //     HighSiteSpinorSet;
    // Same memory layout and on-disk format as the commented spelling above
    // (differs only by an iScalar wrapper); this spelling matches
    // A2Autils::UnpackBinnedVectors, the CPU-view unpack that avoids the
    // peekLorentz device-stack overflow on GPU builds.
    typedef iVector<typename FImpl::SiteSpinor, lowBinSize>  LowSiteSpinorSet;
    typedef iVector<typename FImpl::SiteSpinor, highBinSize> HighSiteSpinorSet;
public:
    // constructor
    TLoadBinnedA2AVecsV(const std::string name);
    // destructor
    virtual ~TLoadBinnedA2AVecsV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};
MODULE_REGISTER_TMP(LoadBinnedA2AVecsV100_96,
                    ARG(TLoadBinnedA2AVecsV<FIMPL, 100, 96>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecsV100_192,
                    ARG(TLoadBinnedA2AVecsV<FIMPL, 100, 192>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecsV200_96,
                    ARG(TLoadBinnedA2AVecsV<FIMPL, 200, 96>), MIO);
MODULE_REGISTER_TMP(LoadBinnedA2AVecsV200_192,
                    ARG(TLoadBinnedA2AVecsV<FIMPL, 200, 192>), MIO);
MODULE_REGISTER_TMP(LoadZBinnedA2AVecsV100_96,
                    ARG(TLoadBinnedA2AVecsV<ZFIMPL, 100, 96>), MIO);
MODULE_REGISTER_TMP(LoadZBinnedA2AVecsV100_192,
                    ARG(TLoadBinnedA2AVecsV<ZFIMPL, 100, 192>), MIO);
MODULE_REGISTER_TMP(LoadZBinnedA2AVecsV200_96,
                    ARG(TLoadBinnedA2AVecsV<ZFIMPL, 200, 96>), MIO);
MODULE_REGISTER_TMP(LoadZBinnedA2AVecsV200_192,
                    ARG(TLoadBinnedA2AVecsV<ZFIMPL, 200, 192>), MIO);
/******************************************************************************
 *                  TLoadBinnedA2AVecsV implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
TLoadBinnedA2AVecsV<FImpl, lowBinSize, highBinSize>::
TLoadBinnedA2AVecsV(const std::string name)
: Module<LoadBinnedA2AVecsVPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
std::vector<std::string>
TLoadBinnedA2AVecsV<FImpl, lowBinSize, highBinSize>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl, int lowBinSize, int highBinSize>
std::vector<std::string>
TLoadBinnedA2AVecsV<FImpl, lowBinSize, highBinSize>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
void TLoadBinnedA2AVecsV<FImpl, lowBinSize, highBinSize>::setup(void)
{
    assert(par().lowSize % lowBinSize == 0);
    assert(par().highSize % highBinSize == 0);
    assert((par().lowSize == 0) or (!par().lowFilestem.empty()));
    assert(!par().highFileStems.empty());
    assert(par().highSize > 0);
    assert(par().nHit > 0);
    // envCreate(std::vector<FermionField>, getName(), 1,
    //           par().lowSize
    //           + par().highSize*par().highFileStems.size(),
    //           envGetGrid(FermionField));
    // Emplace with CpuWrite instead of fill-constructing, so the full
    // fermion field set is not also allocated device-side up front.
    unsigned int total = par().lowSize
                         + par().highSize*par().highFileStems.size();
    auto         grid  = envGetGrid(FermionField);

    envCreate(std::vector<FermionField>, getName(), 1, 0, grid);
    auto &v = envGet(std::vector<FermionField>, getName());
    v.reserve(total);
    for (unsigned int i = 0; i < total; ++i)
    {
        v.emplace_back(grid, CpuWrite);
    }
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int lowBinSize, int highBinSize>
void TLoadBinnedA2AVecsV<FImpl, lowBinSize, highBinSize>::execute(void)
{
    auto &v      = envGet(std::vector<FermionField>, getName());
    int  lowNb   = par().lowSize / lowBinSize;
    int  highNb  = par().highSize / highBinSize;
    Real norm    = 1.0 / par().nHit;

    assert(v.size() ==
           par().lowSize
           + par().highSize*par().highFileStems.size());

    // Low modes
    if (par().lowSize > 0)
    {
        // Lattice<LowSiteSpinorSet> bvec(
        //     envGetGrid(Lattice<LowSiteSpinorSet>));
        // One-element vector so the bin can be handed to
        // A2Autils::UnpackBinnedVectors; still only one bin resident.
        std::vector<Lattice<LowSiteSpinorSet>> bvec(
            1, envGetGrid(Lattice<LowSiteSpinorSet>));
        LOG(Message) << "Loading " << par().lowSize
                     << " low-mode V vectors from "
                     << lowNb << " files of "
                     << lowBinSize << " binned vectors"
                     << std::endl;
        if (par().multiFile)
        {
            for (int i = 0; i < lowNb; ++i)
            {
                A2AVectorsIo::readElement(par().lowFilestem, bvec[0], i,
                                          vm().getTrajectory());
                // for (int j = 0; j < lowBinSize; ++j)
                // {
                //     v[i*lowBinSize + j] = peekLorentz(bvec, j);
                // }
                A2Autils<FImpl>::template UnpackBinnedVectors<lowBinSize>(
                    v, i*lowBinSize, bvec);
            }
        }
        else
        {
            ScidacReader binReader;

            A2AVectorsIo::openReader(binReader, par().lowFilestem,
                                     vm().getTrajectory());
            for (int i = 0; i < lowNb; ++i)
            {
                A2AVectorsIo::readRecord(binReader, bvec[0], i);
                // for (int j = 0; j < lowBinSize; ++j)
                // {
                //     v[i*lowBinSize + j] = peekLorentz(bvec, j);
                // }
                A2Autils<FImpl>::template UnpackBinnedVectors<lowBinSize>(
                    v, i*lowBinSize, bvec);
            }
            binReader.close();
        }
    }

    // High modes
    // Lattice<HighSiteSpinorSet> bvec(
    //     envGetGrid(Lattice<HighSiteSpinorSet>));
    std::vector<Lattice<HighSiteSpinorSet>> bvec(
        1, envGetGrid(Lattice<HighSiteSpinorSet>));
    for (unsigned int is = 0; is < par().highFileStems.size(); ++is)
    {
        unsigned int offset = par().lowSize + is*par().highSize;

        LOG(Message) << "Loading " << par().highSize
                     << " high-mode V vectors from "
                     << par().highFileStems[is] << " in "
                     << highNb << " files of "
                     << highBinSize << " binned vectors"
                     << std::endl;
        if (par().multiFile)
        {
            for (int i = 0; i < highNb; ++i)
            {
                A2AVectorsIo::readElement(par().highFileStems[is], bvec[0], i,
                                          vm().getTrajectory());
                // for (int j = 0; j < highBinSize; ++j)
                // {
                //     v[offset + i*highBinSize + j]
                //         = norm*peekLorentz(bvec, j);
                // }
                A2Autils<FImpl>::template UnpackBinnedVectors<highBinSize>(
                    v, offset + i*highBinSize, bvec);
            }
        }
        else
        {
            ScidacReader binReader;

            A2AVectorsIo::openReader(binReader, par().highFileStems[is],
                                     vm().getTrajectory());
            for (int i = 0; i < highNb; ++i)
            {
                A2AVectorsIo::readRecord(binReader, bvec[0], i);
                // for (int j = 0; j < highBinSize; ++j)
                // {
                //     v[offset + i*highBinSize + j]
                //         = norm*peekLorentz(bvec, j);
                // }
                A2Autils<FImpl>::template UnpackBinnedVectors<highBinSize>(
                    v, offset + i*highBinSize, bvec);
            }
            binReader.close();
        }
        // 1/nHit hit-average normalization, applied to the whole stem's
        // block after the unpack instead of per element above.
        for (unsigned int k = offset; k < offset + par().highSize; ++k)
        {
            v[k] = norm*v[k];
        }
    }
}

END_MODULE_NAMESPACE
END_HADRONS_NAMESPACE
#endif // Hadrons_MIO_LoadBinnedA2AVecsV_hpp_
