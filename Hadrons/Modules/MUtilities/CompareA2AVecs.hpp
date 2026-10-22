/*
 * CompareA2AVecs.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MUtilities_CompareA2AVecs_hpp_
#define Hadrons_MUtilities_CompareA2AVecs_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Exact elementwise comparison of two A2A vector arrays, for regression    *
 *  tests of paths that should move bytes without arithmetic (e.g. the       *
 *  split-grid loader against the full-grid one). Logs how many vectors      *
 *  differ, the largest norm2(left[i] - right[i]) with its index and         *
 *  norm2(left[i]) for scale, and throws if any difference is nonzero --     *
 *  there is no tolerance, since a byte-moving path has no roundoff.         *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class CompareA2AVecsPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CompareA2AVecsPar,
                                    std::string, left,
                                    std::string, right);
};

template <typename FImpl>
class TCompareA2AVecs: public Module<CompareA2AVecsPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TCompareA2AVecs(const std::string name);
    // destructor
    virtual ~TCompareA2AVecs(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CompareA2AVecs, TCompareA2AVecs<FIMPL>, MUtilities);

/******************************************************************************
 *                     TCompareA2AVecs implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TCompareA2AVecs<FImpl>::TCompareA2AVecs(const std::string name)
: Module<CompareA2AVecsPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TCompareA2AVecs<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().left, par().right};

    return in;
}

template <typename FImpl>
std::vector<std::string> TCompareA2AVecs<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCompareA2AVecs<FImpl>::setup(void)
{
    envTmpLat(FermionField, "diff");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TCompareA2AVecs<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    envGetTmp(FermionField, diff);

    if (left.size() != right.size())
    {
        HADRONS_ERROR(Size, "'" + par().left + "' has " + std::to_string(left.size())
                      + " vectors, '" + par().right + "' has "
                      + std::to_string(right.size()));
    }

    unsigned int nDiffer = 0, nZero = 0, worst = 0;
    RealD        maxDiff = 0., minNorm = 0., maxNorm = 0.;

    for (unsigned int i = 0; i < left.size(); ++i)
    {
        diff = left[i] - right[i];
        RealD d = norm2(diff);
        RealD n = norm2(left[i]);
        if (d != 0.)
        {
            ++nDiffer;
        }
        if (d > maxDiff)
        {
            maxDiff = d;
            worst   = i;
        }
        if (n == 0.)
        {
            ++nZero;
        }
        minNorm = (i == 0) ? n : std::min(minNorm, n);
        maxNorm = std::max(maxNorm, n);
    }
    // Zero-norm inputs would make an all-zero difference trivial, hence nZero.
    LOG(Message) << "'" << par().left << "' vs '" << par().right << "': "
                 << nDiffer << " of " << left.size() << " vectors differ" << std::endl;
    LOG(Message) << "max norm2(diff) = " << maxDiff << " at index " << worst
                 << ", norm2(left) = " << norm2(left[worst]) << std::endl;
    LOG(Message) << "norm2(left) in [" << minNorm << ", " << maxNorm << "], "
                 << nZero << " zero-norm vectors" << std::endl;
    if (nDiffer > 0)
    {
        HADRONS_ERROR(Implementation, "A2A vector arrays differ");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_CompareA2AVecs_hpp_
