/*
 * A2AHighModeBatchedV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MSolver_A2AHighModeBatchedV_hpp_
#define Hadrons_MSolver_A2AHighModeBatchedV_hpp_
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
 * Batched counterpart of MSolver::A2AHighModeV. Identical output; the only
 * difference is that batchSize sources are handed to the solver in one call
 * rather than one at a time, through Solver's vector entry point.
 *
 *   batchSize  sources per solver call; 0 and 1 both mean no batching. Must
 *              divide the number of vectors computed by this run.
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
 * bit-identical to A2AHighModeV at any batchSize.
 *
 * Each V vector is written as soon as its batch is solved, one fermion field
 * per file at <output>.<traj>/elem<i>.bin, i the noise vector index. An empty
 * output computes without writing.
 *
 *   ioSplit  split-grid processor layout for writing a batch (e.g. "2 2 2 4";
 *            empty or absent = one file at a time from the full grid). Each
 *            partition writes one vector of the batch to its own file at once,
 *            see A2AVectorsIo::writeBlock. An IO layout only: it is unrelated to
 *            any split the solver uses internally (the srun's
 *            --batched-solver-split), which this module never sees. Timer
 *            "V I/O" is the whole write of each batch.
 *
 * Vector range: initVec and nVec select the noise vectors
 * [initVec, initVec + nVec) computed by this run, so one hit can be split
 * across several jobs (e.g. to fit a walltime cap) or resumed after a failure
 * from the first missing element.
 *
 *   initVec  first noise vector to compute (default 0)
 *   nVec     number of noise vectors to compute; 0 (the default) means all
 *            vectors from initVec to the end of the hit
 *
 * Indices are absolute: vector i is always written as element i, whichever
 * run computes it. Nothing is skipped automatically: which vectors a run
 * computes is exactly what the XML says.
 */
BEGIN_MODULE_NAMESPACE(MSolver)
class A2AHighModeBatchedVPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AHighModeBatchedVPar,
                                  std::string,  noise,
                                  std::string,  action,
                                  std::string,  solver,
                                  std::string,  output,
                                  unsigned int, initVec,
                                  unsigned int, nVec,
                                  unsigned int, batchSize,
                                  std::string,  ioSplit);
};
template <typename FImpl>
class TA2AHighModeBatchedV : public Module<A2AHighModeBatchedVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AHighModes<FImpl> A2A;
public:
    // constructor
    TA2AHighModeBatchedV(const std::string name);
    // destructor
    virtual ~TA2AHighModeBatchedV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    struct VecRange
    {
        unsigned int first;
        unsigned int count;
    };
    // noise vectors computed by this run, validated against the noise size
    VecRange vecRange(unsigned int fermSize);
    // sources per solver call, validated against the vector range
    unsigned int batchSize(const VecRange &range);
};
MODULE_REGISTER_TMP(A2AHighModeBatchedV, ARG(TA2AHighModeBatchedV<FIMPL>), MSolver);
/******************************************************************************
 *                    TA2AHighModeBatchedV implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AHighModeBatchedV<FImpl>::TA2AHighModeBatchedV(const std::string name)
: Module<A2AHighModeBatchedVPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AHighModeBatchedV<FImpl>::getInput(void)
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

template <typename FImpl>
std::vector<std::string> TA2AHighModeBatchedV<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}
// vector range ////////////////////////////////////////////////////////////////
template <typename FImpl>
typename TA2AHighModeBatchedV<FImpl>::VecRange
TA2AHighModeBatchedV<FImpl>::vecRange(unsigned int fermSize)
{
    VecRange range;

    range.first = par().initVec;
    if (range.first >= fermSize)
    {
        HADRONS_ERROR(Argument, "initVec = " + std::to_string(range.first)
                      + " is outside the " + std::to_string(fermSize)
                      + " vectors of noise '" + par().noise + "'");
    }
    range.count = (par().nVec == 0) ? fermSize - range.first : par().nVec;
    if (range.first + range.count > fermSize)
    {
        HADRONS_ERROR(Argument, "vectors [" + std::to_string(range.first) + ", "
                      + std::to_string(range.first + range.count)
                      + ") exceed the " + std::to_string(fermSize)
                      + " vectors of noise '" + par().noise + "'");
    }

    return range;
}
// batch size //////////////////////////////////////////////////////////////////
template <typename FImpl>
unsigned int TA2AHighModeBatchedV<FImpl>::batchSize(const VecRange &range)
{
    unsigned int nBatch = (par().batchSize == 0) ? 1 : par().batchSize;

    if (range.count % nBatch != 0)
    {
        HADRONS_ERROR(Argument, "batchSize = " + std::to_string(nBatch)
                      + " does not divide the " + std::to_string(range.count)
                      + " vectors computed by this run");
    }

    return nBatch;
}
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AHighModeBatchedV<FImpl>::setup(void)
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
    unsigned int nBatch  = batchSize(vecRange(noise.fermSize()));
    envTmp(std::vector<FermionField>, "noiseBatch", 1, nBatch,
           envGetGrid(FermionField));
    envTmp(std::vector<FermionField>, "vBatch", 1, nBatch,
           envGetGrid(FermionField));
    if (Ls > 1)
    {
        envTmp(std::vector<FermionField>, "f5Batch", Ls, nBatch,
               envGetGrid(FermionField, Ls));
    }
    envTmp(A2A, "a2a", 1, action, solver);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AHighModeBatchedV<FImpl>::execute(void)
{
    auto         &noise  = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    int          Ls      = env().getObjectLs(par().action);
    VecRange     range   = vecRange(noise.fermSize());
    unsigned int nBatch  = batchSize(range);
    unsigned int ihBegin = range.first;
    unsigned int ihEnd   = range.first + range.count;
    envGetTmp(std::vector<FermionField>, noiseBatch);
    envGetTmp(std::vector<FermionField>, vBatch);
    envGetTmp(A2A, a2a);
    LOG(Message) << "Computing high-mode part of all-to-all V vectors "
                 << "using noise '" << par().noise << "' ("
                 << noise.fermSize() << " noise vectors)" << std::endl;
    LOG(Message) << "Noise vectors [" << ihBegin << ", " << ihEnd << ")"
                 << std::endl;
    LOG(Message) << "Solving " << nBatch << " source(s) per solver call"
                 << std::endl;
    auto grid = envGetGrid(FermionField);
    std::unique_ptr<GridCartesian> splitGrid;
    if (!par().output.empty() && !par().ioSplit.empty())
    {
        splitGrid.reset(new GridCartesian(grid->FullDimensions(), grid->_simd_layout,
                                          strToVec<int>(par().ioSplit), *grid));
        LOG(Message) << "Split write: "
                     << grid->ProcessorCount()/splitGrid->ProcessorCount()
                     << " partitions of layout " << par().ioSplit << std::endl;
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
        if (!par().output.empty())
        {
            startTimer("V I/O");
            A2AVectorsIo::writeBlock(par().output, vBatch, 0, nBatch,
                                     vm().getTrajectory(), splitGrid.get(), ihBase);
            stopTimer("V I/O");
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MSolver_A2AHighModeBatchedV_hpp_
