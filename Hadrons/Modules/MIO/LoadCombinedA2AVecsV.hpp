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
 *  All files are one fermion field per file, <stem>.<traj>/elem<i>.bin,     *
 *  as written by MSolver::A2AHighModeV and MUtilities::A2ALowModes(Coarse). *
 *  The per-hit high filestem is 'highStem + highExtensions[h]' -- e.g.      *
 *  highStem = ".../vw/" and highExtensions = {"l0_v", "l1_v", ...} -- so    *
 *  the common directory only has to be written once rather than once per   *
 *  hit.                                                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCombinedA2AVecsVPar: Serializable
{
public:
    // nHit: 1/nHit hit-average normalization applied to the high-mode blocks
    // after loading -- the multi-hit estimator averages over hits, and by
    // convention the whole factor lives on the V side while W stays raw
    // unit-modulus noise. Set it to the hit
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

template <typename FImpl>
class TLoadCombinedA2AVecsV: public Module<LoadCombinedA2AVecsVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
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

MODULE_REGISTER_TMP(LoadCombinedA2AVecsV, ARG(TLoadCombinedA2AVecsV<FIMPL>), MIO);

/******************************************************************************
 *                  TLoadCombinedA2AVecsV implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadCombinedA2AVecsV<FImpl>::TLoadCombinedA2AVecsV(const std::string name)
: Module<LoadCombinedA2AVecsVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadCombinedA2AVecsV<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadCombinedA2AVecsV<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadCombinedA2AVecsV<FImpl>::setup(void)
{
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
template <typename FImpl>
void TLoadCombinedA2AVecsV<FImpl>::execute(void)
{
    auto &out = envGet(std::vector<FermionField>, getName());
    unsigned int offset = 0;

    if (par().nLow > 0)
    {
        LOG(Message) << "Loading " << par().nLow << " low-mode A2A vectors from '"
                     << par().lowFilestem << "'" << std::endl;
        for (unsigned int i = 0; i < par().nLow; ++i)
        {
            A2AVectorsIo::readElement(par().lowFilestem, out[offset + i], i,
                                      vm().getTrajectory());
        }
        offset += par().nLow;
    }

    Real norm = (par().nHit > 0) ? (1.0/par().nHit) : 1.0;

    for (unsigned int h = 0; h < par().highExtensions.size(); ++h)
    {
        std::string filestem = par().highStem + par().highExtensions[h];

        LOG(Message) << "Loading " << par().highSize << " high-mode A2A vectors from '"
                     << filestem << "'" << std::endl;
        for (unsigned int i = 0; i < par().highSize; ++i)
        {
            A2AVectorsIo::readElement(filestem, out[offset + i], i,
                                      vm().getTrajectory());
        }
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
