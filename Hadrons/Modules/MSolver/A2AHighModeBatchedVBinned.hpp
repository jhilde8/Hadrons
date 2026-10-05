/*
 * A2AHighModeBatchedVBinned.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn O hOgain <fionnoh@gmail.com>
 * Author: fionnoh <fionnoh@gmail.com>
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
#ifndef Hadrons_MSolver_A2AHighModeBatchedVBinned_hpp_
#define Hadrons_MSolver_A2AHighModeBatchedVBinned_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *            Create high-mode all-to-all V vectors, batched solves           *
 ******************************************************************************/
/*
 * Batched counterpart of MSolver::A2AHighModeVBinned. Identical output; the
 * only difference is that batchSize sources are handed to the solver in one
 * call rather than one at a time, through Solver's vector entry point.
 *
 *   batchSize  sources per solver call; 0 and 1 both mean no batching. Must
 *              divide binSize, so a batch never straddles a bin boundary.
 *
 * What the batch buys depends entirely on the solver it is pointed at. The
 * vector SchurRedBlackBase::operator() builds every red-black source first,
 * then calls its guesser ONCE for the whole batch, so a deflation guesser sees
 * all sources together. The solve itself is only batched if the solver's
 * OperatorFunction overrides the vector operator() -- which
 * MSolver::MixedPrecisionRBPrecCG does not, so with that solver the batch is
 * deflated together and then solved one source at a time by
 * OperatorFunction's base-class loop. That combination is the point of this
 * module: batched deflation with provably unchanged solves, so output must be
 * bit-identical to A2AHighModeVBinned at any batchSize.
 *
 * Bin range: a hit of fermSize noise vectors is written as fermSize/binSize
 * bins. initBin and nBin select the bins [initBin, initBin + nBin) computed by
 * this run, so one hit can be split across several jobs (e.g. to fit a
 * walltime cap) or resumed after a failure from the first incomplete bin.
 *
 *   initBin  first bin to compute (default 0)
 *   nBin     number of bins to compute; 0 (the default) means all bins from
 *            initBin to the end of the hit
 *
 * Bin indices are absolute: bin b is always written as element b, whichever
 * run computes it. Any range other than the full hit requires multiFile=true,
 * since a single-file run over a partial range would leave a truncated file.
 * Nothing is skipped automatically: which bins a run computes is exactly what
 * the XML says.
 */
BEGIN_MODULE_NAMESPACE(MSolver)
class A2AHighModeBatchedVBinnedPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AHighModeBatchedVBinnedPar,
                                  std::string,  noise,
                                  std::string,  action,
                                  std::string,  solver,
                                  std::string,  output,
                                  bool,         multiFile,
                                  unsigned int, initBin,
                                  unsigned int, nBin,
                                  unsigned int, batchSize);
};
template <typename FImpl, int binSize>
class TA2AHighModeBatchedVBinned : public Module<A2AHighModeBatchedVBinnedPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AHighModes<FImpl> A2A;
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2AHighModeBatchedVBinned(const std::string name);
    // destructor
    virtual ~TA2AHighModeBatchedVBinned(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    struct BinRange
    {
        unsigned int first;
        unsigned int count;
    };
    // bins computed by this run, validated against the noise size
    BinRange binRange(unsigned int fermSize);
    // sources per solver call, validated against binSize
    unsigned int batchSize(void);
};
MODULE_REGISTER_TMP(A2AHighModeBatchedVBinned12,
    ARG(TA2AHighModeBatchedVBinned<FIMPL, 12>), MSolver);
MODULE_REGISTER_TMP(A2AHighModeBatchedVBinned16,
    ARG(TA2AHighModeBatchedVBinned<FIMPL, 16>), MSolver);
MODULE_REGISTER_TMP(A2AHighModeBatchedVBinned96,
    ARG(TA2AHighModeBatchedVBinned<FIMPL, 96>), MSolver);
MODULE_REGISTER_TMP(A2AHighModeBatchedVBinned128,
    ARG(TA2AHighModeBatchedVBinned<FIMPL, 128>), MSolver);
MODULE_REGISTER_TMP(A2AHighModeBatchedVBinned192,
    ARG(TA2AHighModeBatchedVBinned<FIMPL, 192>), MSolver);
/******************************************************************************
 *                 TA2AHighModeBatchedVBinned implementation                  *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
TA2AHighModeBatchedVBinned<FImpl, binSize>::TA2AHighModeBatchedVBinned(const std::string name)
: Module<A2AHighModeBatchedVBinnedPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int binSize>
std::vector<std::string> TA2AHighModeBatchedVBinned<FImpl, binSize>::getInput(void)
{
    // MSolver::MixedPrecisionRBPrecCG always creates both a plain and a
    // "_subtract" solver object regardless of whether an outer guesser was
    // configured, so both are safe to list here. setup() below picks
    // whichever one actually has a guesser attached (Solver::hasGuesser()) --
    // subtracting a guess only makes sense, and is only checkerboard-safe,
    // when a real (non-Zero) outer guesser produced that guess.
    std::vector<std::string> in = {par().action, par().solver,
                                   par().solver + "_subtract", par().noise};
    return in;
}

template <typename FImpl, int binSize>
std::vector<std::string> TA2AHighModeBatchedVBinned<FImpl, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}
// bin range ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
typename TA2AHighModeBatchedVBinned<FImpl, binSize>::BinRange
TA2AHighModeBatchedVBinned<FImpl, binSize>::binRange(unsigned int fermSize)
{
    unsigned int nTotal = fermSize/binSize;
    BinRange     range;

    range.first = par().initBin;
    if (range.first >= nTotal)
    {
        HADRONS_ERROR(Argument, "initBin = " + std::to_string(range.first)
                      + " is outside the " + std::to_string(nTotal)
                      + " bins of noise '" + par().noise + "'");
    }
    range.count = (par().nBin == 0) ? nTotal - range.first : par().nBin;
    if (range.first + range.count > nTotal)
    {
        HADRONS_ERROR(Argument, "bins [" + std::to_string(range.first) + ", "
                      + std::to_string(range.first + range.count)
                      + ") exceed the " + std::to_string(nTotal)
                      + " bins of noise '" + par().noise + "'");
    }
    if ((range.count < nTotal) && !par().multiFile)
    {
        HADRONS_ERROR(Argument, "a partial bin range (initBin/nBin) requires "
                      "multiFile = true");
    }

    return range;
}
// batch size //////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
unsigned int TA2AHighModeBatchedVBinned<FImpl, binSize>::batchSize(void)
{
    unsigned int nBatch = (par().batchSize == 0) ? 1 : par().batchSize;

    if (binSize % nBatch != 0)
    {
        HADRONS_ERROR(Argument, "batchSize = " + std::to_string(nBatch)
                      + " does not divide the bin size "
                      + std::to_string(binSize));
    }

    return nBatch;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TA2AHighModeBatchedVBinned<FImpl, binSize>::setup(void)
{
    auto &noise          = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto &action         = envGet(FMat, par().action);
    auto &solverPlain    = envGet(Solver, par().solver);
    auto &solverSubtract = envGet(Solver, par().solver + "_subtract");
    // Only take the guess-subtracting path if a real outer guesser was
    // configured, so the stored high-mode V is the low-mode-deflated remainder
    // and does not double count low modes included exactly from the eigenpack.
    // With no guesser, "_subtract" would subtract a ZeroGuesser guess whose
    // Checkerboard tag is never set to match the RB solution, tripping the
    // checkerboard-consistency assert in Lattice_ET.h.
    auto         &solver = solverSubtract.hasGuesser() ? solverSubtract : solverPlain;
    int          Ls      = env().getObjectLs(par().action);
    unsigned int nBatch  = batchSize();
    assert(noise.fermSize() % binSize == 0);
    binRange(noise.fermSize());
    envTmp(std::vector<FermionField>, "noiseBatch", 1, nBatch,
           envGetGrid(FermionField));
    envTmp(std::vector<FermionField>, "vBatch", 1, nBatch,
           envGetGrid(FermionField));
    envTmpLat(Lattice<SiteSpinorSet>, "vBin");
    if (Ls > 1)
    {
        envTmp(std::vector<FermionField>, "f5Batch", Ls, nBatch,
               envGetGrid(FermionField, Ls));
    }
    envTmp(A2A, "a2a", 1, action, solver);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TA2AHighModeBatchedVBinned<FImpl, binSize>::execute(void)
{
    auto         &noise  = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    int          Ls      = env().getObjectLs(par().action);
    BinRange     range   = binRange(noise.fermSize());
    unsigned int nBatch  = batchSize();
    unsigned int ihBegin = binSize*range.first;
    unsigned int ihEnd   = binSize*(range.first + range.count);
    envGetTmp(std::vector<FermionField>, noiseBatch);
    envGetTmp(std::vector<FermionField>, vBatch);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(A2A, a2a);
    ScidacWriter vWriter(vBin.Grid()->IsBoss());
    LOG(Message) << "Computing high-mode part of all-to-all V vectors "
                 << "using noise '" << par().noise << "' ("
                 << noise.fermSize() << " noise vectors)" << std::endl;
    LOG(Message) << "Bins [" << range.first << ", "
                 << range.first + range.count << ") of "
                 << noise.fermSize()/binSize << " (noise vectors ["
                 << ihBegin << ", " << ihEnd << "))" << std::endl;
    LOG(Message) << "Solving " << nBatch << " source(s) per solver call"
                 << std::endl;
    if ((!par().output.empty()) && (!par().multiFile))
    {
        A2AVectorsIo::openWriter(vWriter, par().output, vBin.Grid(),
                                 vm().getTrajectory());
    }
    // High modes
    for (unsigned int ihBase = ihBegin; ihBase < ihEnd; ihBase += nBatch)
    {
        LOG(Message) << "V vectors [" << ihBase << ", " << ihBase + nBatch
                     << ") (stochastic modes)" << std::endl;
        // getFerm returns a reference to one shared scratch field, so every
        // source has to be copied out before the next call overwrites it.
        startTimer("V noise");
        for (unsigned int j = 0; j < nBatch; ++j)
        {
            noiseBatch[j] = noise.getFerm(ihBase + j);
        }
        stopTimer("V noise");
        startTimer("V high mode");
        if (Ls == 1)
        {
            a2a.makeHighModeV(vBatch, noiseBatch);
        }
        else
        {
            envGetTmp(std::vector<FermionField>, f5Batch);
            a2a.makeHighModeV5D(vBatch, f5Batch, noiseBatch);
        }
        stopTimer("V high mode");
        for (unsigned int j = 0; j < nBatch; ++j)
        {
            pokeLorentz(vBin, vBatch[j], (ihBase + j) % binSize);
        }

        if (((ihBase + nBatch) % binSize == 0) && (!par().output.empty()))
        {
            unsigned int ib = (ihBase + nBatch - 1)/binSize;
            startTimer("V I/O");
            if (par().multiFile)
            {
                A2AVectorsIo::writeElement(par().output, vBin, ib,
                                           vm().getTrajectory());
            }
            else
            {
                A2AVectorsIo::writeRecord(vWriter, vBin, ib);
            }
            stopTimer("V I/O");
        }
    }
    if ((!par().output.empty()) && (!par().multiFile))
    {
        vWriter.close();
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MSolver_A2AHighModeBatchedVBinned_hpp_
