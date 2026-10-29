/*
 * A2ALowModesCoarseBLAS.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
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
#ifndef Hadrons_MUtilities_A2ALowModesCoarseBLAS_hpp_
#define Hadrons_MUtilities_A2ALowModesCoarseBLAS_hpp_

#include <memory>

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  MUtilities::A2ALowModesCoarse with the promotion done by GEMM, one batch *
 *  at a time. Same parameters, same output files.                           *
 *                                                                            *
 *  A2ALowModesCoarse promotes each mode with blockPromote, which streams    *
 *  the whole fine basis once per mode. Here Grid's MultiRHSBlockProject     *
 *  repacks the fine basis once, at the start of execute, into a             *
 *  block-major device buffer (BLAS_V), and then promotes a whole batch of   *
 *  coarse vectors with one batched GEMM over coarse sites: the basis is     *
 *  streamed once per batch. The batch is the write batch, one mode per      *
 *  ioSplit partition; without ioSplit it is one mode, and the GEMM has a    *
 *  single right-hand side.                                                  *
 *                                                                            *
 *  Device memory Hadrons' profiler does not see (deviceVector is not part   *
 *  of the evictable cache, so it counts against what --device-mem leaves):  *
 *                                                                            *
 *    BLAS_V  the fine basis again, ~7.5 GiB per GCD for 200 single          *
 *            precision vectors at 64 nodes; logged at import                *
 *    BLAS_F  one batch of fine vectors, single precision                    *
 *                                                                            *
 *  Both are freed at the end of execute. The batch itself (fBatch, single   *
 *  precision fine vectors) is an envTmp and is seen.                        *
 *                                                                            *
 *  The first GEMM in a process loads rocBLAS's kernels, ~70 s cold at 256   *
 *  ranks; that lands on the first batch's "Promote".                        *
 *                                                                            *
 *  Not bit-compatible with A2ALowModesCoarse: a GEMM and a loop of axpy sum *
 *  in different orders. Compare with a tolerance (Test_low_modes_blas).     *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MUtilities)

class A2ALowModesCoarseBLASPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ALowModesCoarseBLASPar,
                                    std::string, eigenPack,
                                    std::string, action,
                                    std::string, output,
                                    std::string, ioSplit);
};

template <typename FImpl, typename FImplPack, int nBasis>
class TA2ALowModesCoarseBLAS : public Module<A2ALowModesCoarseBLASPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    FERM_TYPE_ALIASES(FImplPack, Pack);
    typedef CoarseFermionEigenPack<FImplPack, nBasis> EPack;
    typedef typename EPack::CoarseField               CoarseField;
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl> A2A;
public:
    // constructor
    TA2ALowModesCoarseBLAS(const std::string name);
    // destructor
    virtual ~TA2ALowModesCoarseBLAS(void) {};
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

MODULE_REGISTER_TMP(A2ALowModesCoarseBLAS200,
    ARG(TA2ALowModesCoarseBLAS<FIMPL, FIMPLF, 200>), MUtilities);
MODULE_REGISTER_TMP(ZA2ALowModesCoarseBLAS200,
    ARG(TA2ALowModesCoarseBLAS<ZFIMPL, ZFIMPLF, 200>), MUtilities);

/******************************************************************************
 *                   TA2ALowModesCoarseBLAS implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
TA2ALowModesCoarseBLAS<FImpl, FImplPack, nBasis>::TA2ALowModesCoarseBLAS(const std::string name)
: Module<A2ALowModesCoarseBLASPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
std::vector<std::string> TA2ALowModesCoarseBLAS<FImpl, FImplPack, nBasis>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack, par().action};

    return in;
}

template <typename FImpl, typename FImplPack, int nBasis>
std::vector<std::string> TA2ALowModesCoarseBLAS<FImpl, FImplPack, nBasis>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
void TA2ALowModesCoarseBLAS<FImpl, FImplPack, nBasis>::setup(void)
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
    envTmp(std::vector<FermionFieldPack>, "fBatch", Ls, nBatch_, epack.evec[0].Grid());
    envTmp(std::vector<CoarseField>, "cBatch", Ls, nBatch_, epack.evecCoarse[0].Grid());
    envTmp(FermionField, "evecD", Ls, action.FermionRedBlackGrid());
    envTmp(std::vector<FermionField>, "vBatch", 1, nBatch_, envGetGrid(FermionField));
    envTmp(std::vector<FermionField>, "wBatch", 1, nBatch_, envGetGrid(FermionField));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, typename FImplPack, int nBasis>
void TA2ALowModesCoarseBLAS<FImpl, FImplPack, nBasis>::execute(void)
{
    auto        &epack  = envGet(EPack, par().eigenPack);
    int         Ls      = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);
    envGetTmp(std::vector<FermionFieldPack>, fBatch);
    envGetTmp(std::vector<CoarseField>, cBatch);
    envGetTmp(FermionField, evecD);
    envGetTmp(std::vector<FermionField>, vBatch);
    envGetTmp(std::vector<FermionField>, wBatch);

    LOG(Message) << "Streaming low-mode all-to-all vectors, promoting by GEMM "
                 << "from compressed coarse eigenpack '" << par().eigenPack
                 << "' (" << Nl_ << " low modes, batches of " << nBatch_ << ")"
                 << std::endl;

    MultiRHSBlockProject<FermionFieldPack> blockProj;

    startTimer("Import basis");
    blockProj.Allocate(epack.evec.size(), epack.evec[0].Grid(), epack.evecCoarse[0].Grid());
    blockProj.ImportBasis(epack.evec);
    stopTimer("Import basis");
    LOG(Message) << "Packed " << epack.evec.size() << " fine basis vectors, "
                 << blockProj.DeviceBytes()/1024./1024./1024.
                 << " GiB resident device memory per rank" << std::endl;

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

        // Only the last batch can be short; blockPromote takes nrhs from the size.
        if (n < fBatch.size())
        {
            fBatch.erase(fBatch.begin() + n, fBatch.end());
            cBatch.erase(cBatch.begin() + n, cBatch.end());
        }

        startTimer("Promote");
        for (unsigned int j = 0; j < n; j++)
        {
            cBatch[j] = epack.evecCoarse[il0 + j];
        }
        blockProj.blockPromote(fBatch, cBatch);
        stopTimer("Promote");

        for (unsigned int j = 0; j < n; j++)
        {
            unsigned int il = il0 + j;

            // The GEMM never sets a checkerboard; tag it Odd like epack.evec.
            startTimer("Promote");
            fBatch[j].Checkerboard() = epack.evec[0].Checkerboard();
            precisionChange(evecD, fBatch[j]);
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

#endif // Hadrons_MUtilities_A2ALowModesCoarseBLAS_hpp_
