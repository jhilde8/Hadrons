/*
 * VectorPackRefSlice.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_VectorPackRefSlice_hpp_
#define Hadrons_MUtilities_VectorPackRefSlice_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Pack pointers to one or more offset/count segments of a single existing  *
 *  field array into a std::vector<const Field *>, concatenated in the       *
 *  order the segments are listed.                                          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class VectorPackRefSlicePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(VectorPackRefSlicePar,
                                    std::string,                source,
                                    std::vector<unsigned int>, offsets,
                                    std::vector<unsigned int>, counts);
};

template <typename Field>
class TVectorPackRefSlice: public Module<VectorPackRefSlicePar>
{
public:
    // constructor
    TVectorPackRefSlice(const std::string name);
    // destructor
    virtual ~TVectorPackRefSlice(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(FermionVectorPackRefSlice, TVectorPackRefSlice<FIMPL::FermionField>, MUtilities);

/******************************************************************************
 *                 TVectorPackRefSlice implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Field>
TVectorPackRefSlice<Field>::TVectorPackRefSlice(const std::string name)
: Module<VectorPackRefSlicePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Field>
std::vector<std::string> TVectorPackRefSlice<Field>::getInput(void)
{
    std::vector<std::string> in = {par().source};

    return in;
}

template <typename Field>
std::vector<std::string> TVectorPackRefSlice<Field>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename Field>
DependencyMap TVectorPackRefSlice<Field>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().source, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorPackRefSlice<Field>::setup(void)
{
    if (par().offsets.size() != par().counts.size())
    {
        HADRONS_ERROR(Size, "offsets and counts must have the same length");
    }

    unsigned int total = 0;
    for (auto c: par().counts)
    {
        total += c;
    }

    envCreate(std::vector<const Field *>, getName(), 1, total, nullptr);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Field>
void TVectorPackRefSlice<Field>::execute(void)
{
    auto &src = envGet(std::vector<Field>, par().source);
    auto &vec = envGet(std::vector<const Field *>, getName());

    LOG(Message) << "Packing pointers from " << par().offsets.size()
                 << " segment(s) of '" << par().source << "'" << std::endl;

    unsigned int cursor = 0;
    for (unsigned int i = 0; i < par().offsets.size(); ++i)
    {
        for (unsigned int j = 0; j < par().counts[i]; ++j)
        {
            vec[cursor++] = &src[par().offsets[i] + j];
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_VectorPackRefSlice_hpp_
