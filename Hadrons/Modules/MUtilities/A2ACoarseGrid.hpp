/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MUtilities/A2ACoarseGrid.hpp

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
#ifndef Hadrons_MUtilities_A2ACoarseGrid_hpp_
#define Hadrons_MUtilities_A2ACoarseGrid_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                    Module to load all-to-all vectors                       *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ACoarseGridPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACoarseGridPar,
				    std::string,  fine,
				    std::vector<int>, blockSize,
				    std::vector<int>, offsets,
				    std::string, output);
};

template <typename FImpl>
class TA2ACoarseGrid: public Module<A2ACoarseGridPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TA2ACoarseGrid(const std::string name);
    // destructor
    virtual ~TA2ACoarseGrid(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(A2ACoarseGrid, TA2ACoarseGrid<FIMPL>, MUtilities);

/******************************************************************************
 *                       TA2ACoarseGrid implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2ACoarseGrid<FImpl>::TA2ACoarseGrid(const std::string name)
: Module<A2ACoarseGridPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2ACoarseGrid<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().fine};
    
  return in;
}

template <typename FImpl>
std::vector<std::string> TA2ACoarseGrid<FImpl>::getOutput(void)
{
  std::vector<std::string> out = {getName()};
  
  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ACoarseGrid<FImpl>::setup(void)
{
  for( int mu = 0 ; mu < 4 ; ++mu )
    assert( par().blockSize[mu] > par().offsets[mu]);
  auto &fine   = envGet(std::vector<FermionField>, par().fine);
  envCreate(std::vector<FermionField>, getName(), 1, fine.size(),
              envGetCoarseGrid(FermionField,par().blockSize));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2ACoarseGrid<FImpl>::execute(void)
{
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef iSpinColourVector<scalar_type> SpinColourVector_s;

  auto &fine   = envGet(std::vector<FermionField>, par().fine);
  auto &coarse = envGet(std::vector<FermionField>, getName());

  GridBase *cgrid = coarse[0].Grid();
  for ( int i = 0 ; i < fine.size() ; ++i ) {
    thread_for(idx,cgrid->lSites(),{
      Coordinate Site(4); Coordinate CoarseSite(4);
      cgrid->LocalIndexToLocalCoor(idx,CoarseSite);
      Site[0] = CoarseSite[0] * par().blockSize[0] + par().offsets[0];
      Site[1] = CoarseSite[1] * par().blockSize[1] + par().offsets[1];
      Site[2] = CoarseSite[2] * par().blockSize[2] + par().offsets[2];
      Site[3] = CoarseSite[3] * par().blockSize[3] + par().offsets[3];
      SpinColourVector_s scvec;
      peekLocalSite(scvec, fine[i], Site);
      pokeLocalSite(scvec, coarse[i],CoarseSite);
    });
  }

  // I/O if necessary
  if (!par().output.empty())
  {
    A2AVectorsIo::write(par().output, coarse, true, vm().getTrajectory());
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ACoarseGrid_hpp_
