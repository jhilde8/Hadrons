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
 *  Sparsen an A2A vector array onto the coarse grid of blockSize, keeping   *
 *  the site at offsets in each block, and optionally write it to output     *
 *  (<output>.<traj>/elem<i>.bin).                                           *
 *                                                                            *
 *  modeOffset and modeSize select fine[modeOffset, modeOffset + modeSize)   *
 *  (modeSize 0 = to the end), so one block of a combined array -- the low  *
 *  modes, or one hit -- is sparsened and written as its own stem. Every     *
 *  coarse field is multiplied by scale (required; 1 = unchanged), e.g. nHit *
 *  to write a V block raw after the loader applied 1/nHit.                  *
 *                                                                            *
 *  ioSplit (e.g. "2 2 2 1"; empty or absent = off) writes on split-grid     *
 *  partitions of the coarse grid with that processor layout. A sparsened    *
 *  field is too small to write efficiently from the full grid -- ~25 KB per *
 *  rank on 64I at 512 nodes against the ~4 MB a Lustre write needs -- so    *
 *  Grid_split gathers a batch of fields, one per partition, and each        *
 *  partition writes its own file at once: 12.6 MB per rank on 1-node        *
 *  partitions. Timer "IO" is the whole write.                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ACoarseGridPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ACoarseGridPar,
				    std::string,  fine,
				    unsigned int, modeOffset,
				    unsigned int, modeSize,
				    double, scale,
				    std::vector<int>, blockSize,
				    std::vector<int>, offsets,
				    std::string, output,
				    std::string, ioSplit);
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
  unsigned int size = (par().modeSize == 0) ? fine.size() - par().modeOffset
                                            : par().modeSize;
  if (par().modeOffset + size > fine.size())
  {
    HADRONS_ERROR(Size, "modes [" + std::to_string(par().modeOffset) + ", "
                  + std::to_string(par().modeOffset + size) + ") exceed '"
                  + par().fine + "' (" + std::to_string(fine.size()) + ")");
  }
  envCreate(std::vector<FermionField>, getName(), 1, size,
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
  unsigned int first = par().modeOffset;
  startTimer("Sparsen");
  for ( int i = 0 ; i < coarse.size() ; ++i ) {
    thread_for(idx,cgrid->lSites(),{
      Coordinate Site(4); Coordinate CoarseSite(4);
      cgrid->LocalIndexToLocalCoor(idx,CoarseSite);
      Site[0] = CoarseSite[0] * par().blockSize[0] + par().offsets[0];
      Site[1] = CoarseSite[1] * par().blockSize[1] + par().offsets[1];
      Site[2] = CoarseSite[2] * par().blockSize[2] + par().offsets[2];
      Site[3] = CoarseSite[3] * par().blockSize[3] + par().offsets[3];
      SpinColourVector_s scvec;
      peekLocalSite(scvec, fine[first + i], Site);
      pokeLocalSite(scvec, coarse[i],CoarseSite);
    });
  }
  if (par().scale != 1.0)
  {
    LOG(Message) << "Scaling by " << par().scale << std::endl;
    for (auto &c: coarse)
    {
      c = par().scale*c;
    }
  }
  stopTimer("Sparsen");

  // I/O if necessary
  if (!par().output.empty())
  {
    auto cgridIo = envGetCoarseGrid(FermionField, par().blockSize);

    std::unique_ptr<GridCartesian> splitGrid;
    if (!par().ioSplit.empty())
    {
      splitGrid.reset(new GridCartesian(cgridIo->FullDimensions(), cgridIo->_simd_layout,
                                        strToVec<int>(par().ioSplit), *cgridIo));
      LOG(Message) << "Split write: "
                   << cgridIo->ProcessorCount()/splitGrid->ProcessorCount()
                   << " partitions of layout " << par().ioSplit << std::endl;
    }
    startTimer("IO");
    A2AVectorsIo::writeBlock(par().output, coarse, 0, coarse.size(),
                             vm().getTrajectory(), splitGrid.get());
    stopTimer("IO");
  }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ACoarseGrid_hpp_
