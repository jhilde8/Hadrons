/*
 * A2AHighModeV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Fionn Ó hÓgáin <fionnoh@gmail.com>
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
#ifndef Hadrons_MSolver_A2AHighModeV_hpp_
#define Hadrons_MSolver_A2AHighModeV_hpp_
#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *                    Create high-mode all-to-all V vectors                   *
 ******************************************************************************/
/*
 * Each V vector is written as soon as it is solved, one fermion field per
 * file at <output>.<traj>/elem<i>.bin, i the noise vector index. An empty
 * output computes without writing.
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
class A2AHighModeVPar: Serializable
{
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(A2AHighModeVPar,
                                  std::string,  noise,
                                  std::string,  action,
                                  std::string,  solver,
                                  std::string,  output,
                                  unsigned int, initVec,
                                  unsigned int, nVec);
};
template <typename FImpl>
class TA2AHighModeV : public Module<A2AHighModeVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef A2AHighModes<FImpl> A2A;
public:
    // constructor
    TA2AHighModeV(const std::string name);
    // destructor
    virtual ~TA2AHighModeV(void) {};
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
};
MODULE_REGISTER_TMP(A2AHighModeV,  ARG(TA2AHighModeV<FIMPL>),  MSolver);
MODULE_REGISTER_TMP(ZA2AHighModeV, ARG(TA2AHighModeV<ZFIMPL>), MSolver);
/******************************************************************************
 *                       TA2AHighModeV implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AHighModeV<FImpl>::TA2AHighModeV(const std::string name)
: Module<A2AHighModeVPar>(name)
{}
// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TA2AHighModeV<FImpl>::getInput(void)
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
std::vector<std::string> TA2AHighModeV<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    return out;
}
// vector range ////////////////////////////////////////////////////////////////
template <typename FImpl>
typename TA2AHighModeV<FImpl>::VecRange
TA2AHighModeV<FImpl>::vecRange(unsigned int fermSize)
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
// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AHighModeV<FImpl>::setup(void)
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
    auto &solver = solverSubtract.hasGuesser() ? solverSubtract : solverPlain;
    int  Ls      = env().getObjectLs(par().action);
    vecRange(noise.fermSize());
    envTmpLat(FermionField, "v");
    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
}
// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AHighModeV<FImpl>::execute(void)
{
    auto         &noise  = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    int          Ls      = env().getObjectLs(par().action);
    VecRange     range   = vecRange(noise.fermSize());
    unsigned int ihBegin = range.first;
    unsigned int ihEnd   = range.first + range.count;
    envGetTmp(FermionField, v);
    envGetTmp(A2A, a2a);
    LOG(Message) << "Computing high-mode part of all-to-all V vectors "
                 << "using noise '" << par().noise << "' ("
                 << noise.fermSize() << " noise vectors)" << std::endl;
    LOG(Message) << "Noise vectors [" << ihBegin << ", " << ihEnd << ")"
                 << std::endl;
    // High modes
    for (unsigned int ih = ihBegin; ih < ihEnd; ih++)
    {
        startTimer("V high mode");
        LOG(Message) << "V vector i = " << ih
                     << " (stochastic mode)" << std::endl;
        if (Ls == 1)
        {
            a2a.makeHighModeV(v, noise.getFerm(ih));
        }
        else
        {
            envGetTmp(FermionField, f5);
            a2a.makeHighModeV5D(v, f5, noise.getFerm(ih));
        }
        stopTimer("V high mode");
        if (!par().output.empty())
        {
            startTimer("V I/O");
            A2AVectorsIo::writeElement(par().output, v, ih,
                                       vm().getTrajectory());
            stopTimer("V I/O");
        }
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE
#endif // Hadrons_MSolver_A2AHighModeV_hpp_
