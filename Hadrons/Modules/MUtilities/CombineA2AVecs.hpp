/*
 * CombineA2AVecs.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
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
#ifndef Hadrons_MUtilities_CombineA2AVecs_hpp_
#define Hadrons_MUtilities_CombineA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Splice a low-mode A2A vector array and a high-mode array, both already   *
 *  resident in the environment, into one contiguous array: the first        *
 *  'nLow' elements of 'low' at indices [0, nLow), followed by all of        *
 *  'high' at [nLow, nLow + Nh). This matches the index convention           *
 *  MSolver::A2AVectorsCoarse uses and what MContraction::A2AMesonField      *
 *  expects as 'left'/'right'. This module does no I/O.                      *
 *                                                                            *
 *  Hadrons does not allow two module instances to both produce an object    *
 *  under the same name, so a single low-mode array cannot be mutated in     *
 *  place across hits. Instead, chain the calls: hit 1 combines the          *
 *  pristine low-mode load (e.g. from MIO::LoadBinnedA2AVecs) with hit 1's   *
 *  high modes; hit 2's 'low' input then points at *hit 1's output*, with    *
 *  'nLow' unchanged (its own high-mode tail is simply not copied), and so   *
 *  on. Since each combine call and the previous hit's meson-field module    *
 *  only ever read the previous combined array, never write it, their       *
 *  order relative to each other does not matter, and Hadrons frees each     *
 *  combined array once both have run. At no point does more than one       *
 *  array hold the full low-mode data.                                      *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class CombineA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CombineA2AVecsPar,
                                    std::string,  low,
                                    unsigned int, nLow,
                                    std::string,  high);
};

template <typename FImpl>
class TCombineA2AVecs: public Module<CombineA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCombineA2AVecs(const std::string name);
    // destructor
    virtual ~TCombineA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CombineA2AVecs, TCombineA2AVecs<FIMPL>, MUtilities);

/******************************************************************************
 *                      TCombineA2AVecs implementation                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCombineA2AVecs<FImpl>::TCombineA2AVecs(const std::string name)
: Module<CombineA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCombineA2AVecs<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().low, par().high};

    return in;
}

template <typename FImpl>
std::vector<std::string> TCombineA2AVecs<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCombineA2AVecs<FImpl>::setup(void)
{
    auto &low  = envGet(std::vector<FermionField>, par().low);
    auto &high = envGet(std::vector<FermionField>, par().high);

    if (par().nLow > low.size())
    {
        HADRONS_ERROR(Size, "nLow (" + std::to_string(par().nLow)
                            + ") exceeds the size of '" + par().low + "' ("
                            + std::to_string(low.size()) + ")");
    }

    envCreate(std::vector<FermionField>, getName(), 1, par().nLow + high.size(),
             envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCombineA2AVecs<FImpl>::execute(void)
{
    auto &low  = envGet(std::vector<FermionField>, par().low);
    auto &high = envGet(std::vector<FermionField>, par().high);
    auto &out  = envGet(std::vector<FermionField>, getName());

    unsigned int Nl = par().nLow;
    unsigned int Nh = high.size();

    LOG(Message) << "Combining " << Nl << " low mode(s) from '" << par().low
                 << "' and " << Nh << " high mode(s) from '" << par().high
                 << "' into '" << getName() << "'" << std::endl;

    for (unsigned int i = 0; i < Nl; i++)
    {
        out[i] = low[i];
    }
    for (unsigned int i = 0; i < Nh; i++)
    {
        out[Nl + i] = high[i];
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_CombineA2AVecs_hpp_
