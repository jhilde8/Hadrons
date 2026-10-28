/*
 * A2ALowModesCoarse.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MUtilities_A2ALowModesCoarse_hpp_
#define Hadrons_MUtilities_A2ALowModesCoarse_hpp_

#include <memory>

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Stream low-mode A2A V & W vectors straight to disk directly from a       *
 *  compressed (LC) coarse eigenpack, one fermion field per file at          *
 *  <output>_v.<traj>/elem<i>.bin and <output>_w.<traj>/elem<i>.bin. Low     *
 *  modes are decompressed/promoted on the fly, one batch at a time, so the  *
 *  full Nl-sized V/W arrays are never resident -- only one batch of V and W.*
 *  Coarse-deflation companion to MUtilities::A2ALowModes, which does the    *
 *  same streaming from an uncompressed (exact) eigenpack. Lives in          *
 *  MUtilities rather than MSolver because no solve ever happens here: low   *
 *  modes are built algebraically from the eigenvectors, so unlike           *
 *  MSolver::A2AHighModeV there is no solver dependency at all.              *
 *                                                                            *
 *  ioSplit (e.g. "2 1 2 2"; empty = off) writes each batch on split-grid    *
 *  partitions of that layout, one mode per partition, through               *
 *  A2AVectorsIo::writeBlock. A batch is one mode per partition, so the      *
 *  batch size is the partition count; without ioSplit it is one mode,       *
 *  written from the full grid.                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ALowModesCoarsePar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALowModesCoarsePar,
                                    std::string, eigenPack,
                                    std::string, action,
                                    std::string, output,
                                    std::string, ioSplit);
};

template <typename FImpl, typename FImplPack, int nBasis>
class TA2ALowModesCoarse : public Module<A2ALowModesCoarsePar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    FERM_TYPE_ALIASES(FImplPack, Pack);
    typedef CoarseFermionEigenPack<FImplPack, nBasis> EPack;
    // Concrete class fixed at compile time by HADRONS_DEFAULT_SCHUR (one build
    // per ensemble), so nothing here names a convention.
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2ALowModesCoarse(const std::string name);
    // destructor
    virtual ~TA2ALowModesCoarse(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nl_{0};
    unsigned int nBatch_{1};
};

MODULE_REGISTER_TMP(A2ALowModesCoarse200,
    ARG(TA2ALowModesCoarse<FIMPL, FIMPLF, 200>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModesCoarse200,
    ARG(TA2ALowModesCoarse<ZFIMPL, ZFIMPLF, 200>), MUtilities);

/******************************************************************************
 *                     TA2ALowModesCoarse implementation                      *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
TA2ALowModesCoarse<FImpl, FImplPack, nBasis>::TA2ALowModesCoarse(const std::string name)
: Module<A2ALowModesCoarsePar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
std::vector<std::string> TA2ALowModesCoarse<FImpl, FImplPack, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack, par().action};

    return in;
}

template <typename FImpl, typename FImplPack, int nBasis>
std::vector<std::string> TA2ALowModesCoarse<FImpl, FImplPack, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
void TA2ALowModesCoarse<FImpl, FImplPack, nBasis>::setup(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    auto        &action = envGet(FMat, par().action);
    int         Ls      = env().getObjectLs(par().action);

    if (env().getObjectLs(par().eigenPack) != Ls)
    {
        HADRONS_ERROR(Size, "eigenPack and action Ls mismatch");
    }

    envTmp(A2A, "a2a", 1, action);

    Nl_     = epack.evecCoarse.size();
    nBatch_ = 1;
    if (!par().ioSplit.empty())
    {
        int splitRanks = 1;

        for (int d: strToVec<int>(par().ioSplit))
        {
            splitRanks *= d;
        }
        nBatch_ = envGetGrid(FermionField)->ProcessorCount()/splitRanks;
    }

    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(FermionFieldPack, "evecF", Ls, epack.evec[0].Grid());
    envTmp(FermionField, "evecD", Ls, action.FermionRedBlackGrid());
    envTmp(std::vector<FermionField>, "vBatch", 1, nBatch_, envGetGrid(FermionField));
    envTmp(std::vector<FermionField>, "wBatch", 1, nBatch_, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
void TA2ALowModesCoarse<FImpl, FImplPack, nBasis>::execute(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    int         Ls      = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);
    envGetTmp(FermionFieldPack, evecF);
    envGetTmp(FermionField, evecD);
    envGetTmp(std::vector<FermionField>, vBatch);
    envGetTmp(std::vector<FermionField>, wBatch);

    // blockPromote() never sets its output's checkerboard; tag it Odd like epack.evec.
    evecF.Checkerboard() = epack.evec[0].Checkerboard();

    LOG(Message) << "Streaming low-mode all-to-all vectors, decompressing/"
                 << "promoting on the fly from compressed coarse eigenpack '"
                 << par().eigenPack << "' (" << Nl_ << " low modes)" << std::endl;

    auto grid = envGetGrid(FermionField);
    std::unique_ptr<GridCartesian> splitGrid;
    if (!par().ioSplit.empty())
    {
        splitGrid.reset(new GridCartesian(grid->FullDimensions(), grid->_simd_layout,
                                          strToVec<int>(par().ioSplit), *grid));
        LOG(Message) << "Split write: " << nBatch_ << " partitions of layout "
                     << par().ioSplit << std::endl;
    }

    for (unsigned int il0 = 0; il0 < Nl_; il0 += nBatch_)
    {
        unsigned int n = std::min(nBatch_, Nl_ - il0);

        for (unsigned int j = 0; j < n; j++)
        {
            unsigned int il = il0 + j;

            startTimer("Promote");
            blockPromote(epack.evecCoarse[il], evecF, epack.evec);
            precisionChange(evecD, evecF);
            stopTimer("Promote");

            startTimer("V low mode");
            LOG(Message) << "V vector i = " << il << " (low mode)" << std::endl;
            if (Ls == 1)
            {
                a2a.makeLowModeV(vBatch[j], evecD, epack.evalCoarse[il]);
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a.makeLowModeV5D(vBatch[j], f5, evecD, epack.evalCoarse[il]);
            }
            stopTimer("V low mode");

            startTimer("W low mode");
            LOG(Message) << "W vector i = " << il << " (low mode)" << std::endl;
            if (Ls == 1)
            {
                a2a.makeLowModeW(wBatch[j], evecD, epack.evalCoarse[il]);
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a.makeLowModeW5D(wBatch[j], f5, evecD, epack.evalCoarse[il]);
            }
            stopTimer("W low mode");
        }

        startTimer("V I/O");
        A2AVectorsIo::writeBlock(par().output + "_v", vBatch, 0, n,
                                 vm().getTrajectory(), splitGrid.get(), il0, this);
        stopTimer("V I/O");
        startTimer("W I/O");
        A2AVectorsIo::writeBlock(par().output + "_w", wBatch, 0, n,
                                 vm().getTrajectory(), splitGrid.get(), il0, this);
        stopTimer("W I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MUtilities_A2ALowModesCoarse_hpp_
