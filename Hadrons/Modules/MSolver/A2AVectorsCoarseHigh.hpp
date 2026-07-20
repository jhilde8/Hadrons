/*
 * A2AVectorsCoarseHigh.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MSolver_A2AVectorsCoarseHigh_hpp_
#define Hadrons_MSolver_A2AVectorsCoarseHigh_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Solver.hpp>
#include <Hadrons/A2AVectors.hpp>
#include <Hadrons/DilutedNoise.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  Stream high-mode (stochastic) A2A V & W vectors straight to disk,        *
 *  binned, one bin at a time -- companion to MSolver::A2AVectorsCoarseLow.  *
 *  binSize is the time extent of the ensemble being run (noise is diluted   *
 *  in spin/colour/time, so an NT-sized bin always divides the stochastic    *
 *  vector count evenly). Only one bin's worth of fields is ever resident,   *
 *  instead of holding the full high-mode array (MSolver::A2AVectorsCoarse)  *
 *  and then a second, packed copy of it (MIO::SaveBinnedA2AVecs).           *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSolver)

class A2AVectorsCoarseHighPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AVectorsCoarseHighPar,
                                    std::string, noise,
                                    std::string, action,
                                    std::string, solver,
                                    std::string, output);
};

template <typename FImpl, int binSize>
class TA2AVectorsCoarseHigh : public Module<A2AVectorsCoarseHighPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    typedef HADRONS_DEFAULT_SCHUR_A2A<FImpl>        A2A;
    typedef typename FImpl::SiteSpinor::vector_type vector_type;
    typedef iVector<iVector<iVector<vector_type, Nc>, Ns>, binSize> SiteSpinorSet;
public:
    // constructor
    TA2AVectorsCoarseHigh(const std::string name);
    // destructor
    virtual ~TA2AVectorsCoarseHigh(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Nh_{0};
    unsigned int Nb_{0};
};

MODULE_REGISTER_TMP(A2AVectorsCoarseHigh48,   ARG(TA2AVectorsCoarseHigh<FIMPL, 48>),   MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseHigh64,   ARG(TA2AVectorsCoarseHigh<FIMPL, 64>),   MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseHigh96,   ARG(TA2AVectorsCoarseHigh<FIMPL, 96>),   MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseHigh128,  ARG(TA2AVectorsCoarseHigh<FIMPL, 128>),  MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseHigh144,  ARG(TA2AVectorsCoarseHigh<FIMPL, 144>),  MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseHigh192,  ARG(TA2AVectorsCoarseHigh<FIMPL, 192>),  MSolver);
MODULE_REGISTER_TMP(A2AVectorsCoarseHigh256,  ARG(TA2AVectorsCoarseHigh<FIMPL, 256>),  MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh48,  ARG(TA2AVectorsCoarseHigh<ZFIMPL, 48>),  MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh64,  ARG(TA2AVectorsCoarseHigh<ZFIMPL, 64>),  MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh96,  ARG(TA2AVectorsCoarseHigh<ZFIMPL, 96>),  MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh128, ARG(TA2AVectorsCoarseHigh<ZFIMPL, 128>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh144, ARG(TA2AVectorsCoarseHigh<ZFIMPL, 144>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh192, ARG(TA2AVectorsCoarseHigh<ZFIMPL, 192>), MSolver);
MODULE_REGISTER_TMP(ZA2AVectorsCoarseHigh256, ARG(TA2AVectorsCoarseHigh<ZFIMPL, 256>), MSolver);

/******************************************************************************
 *                   TA2AVectorsCoarseHigh implementation                     *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
TA2AVectorsCoarseHigh<FImpl, binSize>::TA2AVectorsCoarseHigh(const std::string name)
: Module<A2AVectorsCoarseHighPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl, int binSize>
std::vector<std::string> TA2AVectorsCoarseHigh<FImpl, binSize>::getInput(void)
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
std::vector<std::string> TA2AVectorsCoarseHigh<FImpl, binSize>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TA2AVectorsCoarseHigh<FImpl, binSize>::setup(void)
{
    auto        &noise          = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    auto        &action         = envGet(FMat, par().action);
    auto        &solverPlain    = envGet(Solver, par().solver);
    auto        &solverSubtract = envGet(Solver, par().solver + "_subtract");
    // Only take the guess-subtracting path if a real outer guesser was
    // configured. With no guesser, "_subtract" would subtract a ZeroGuesser
    // guess whose Checkerboard tag is never set to match the RB solution,
    // tripping the checkerboard-consistency assert in Lattice_ET.h.
    auto        &solver = solverSubtract.hasGuesser() ? solverSubtract : solverPlain;
    int         Ls      = env().getObjectLs(par().action);

    Nh_ = noise.fermSize();
    if (Nh_ % binSize != 0)
    {
        HADRONS_ERROR(Size, "number of high modes (" + std::to_string(Nh_)
                            + ") is not a multiple of the bin size ("
                            + std::to_string(binSize) + ")");
    }
    Nb_ = Nh_ / binSize;

    if (Ls > 1)
    {
        envTmpLat(FermionField, "f5", Ls);
    }
    envTmp(A2A, "a2a", 1, action, solver);
    envTmp(FermionField, "vTmp", 1, envGetGrid(FermionField));
    envTmp(FermionField, "wTmp", 1, envGetGrid(FermionField));
    envTmp(Lattice<SiteSpinorSet>, "vBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
    envTmp(Lattice<SiteSpinorSet>, "wBin", 1, envGetGrid(Lattice<SiteSpinorSet>));
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl, int binSize>
void TA2AVectorsCoarseHigh<FImpl, binSize>::execute(void)
{
    auto        &noise = envGet(SpinColorDiagonalNoise<FImpl>, par().noise);
    int         Ls     = env().getObjectLs(par().action);

    envGetTmp(A2A, a2a);
    envGetTmp(FermionField, vTmp);
    envGetTmp(FermionField, wTmp);
    envGetTmp(Lattice<SiteSpinorSet>, vBin);
    envGetTmp(Lattice<SiteSpinorSet>, wBin);

    LOG(Message) << "Streaming high-mode all-to-all vectors from noise '"
                 << par().noise << "' (" << Nh_ << " stochastic modes, "
                 << Nb_ << " bins of " << binSize << ")" << std::endl;

    for (unsigned int b = 0; b < Nb_; b++)
    {
        for (unsigned int j = 0; j < binSize; j++)
        {
            unsigned int ih = b * binSize + j;

            startTimer("V high mode");
            LOG(Message) << "V vector i = " << ih << " (high mode)" << std::endl;
            if (Ls == 1)
            {
                a2a.makeHighModeV(vTmp, noise.getFerm(ih));
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a.makeHighModeV5D(vTmp, f5, noise.getFerm(ih));
            }
            pokeLorentz(vBin, vTmp, j);
            stopTimer("V high mode");

            startTimer("W high mode");
            LOG(Message) << "W vector i = " << ih << " (high mode)" << std::endl;
            if (Ls == 1)
            {
                a2a.makeHighModeW(wTmp, noise.getFerm(ih));
            }
            else
            {
                envGetTmp(FermionField, f5);
                a2a.makeHighModeW5D(wTmp, f5, noise.getFerm(ih));
            }
            pokeLorentz(wBin, wTmp, j);
            stopTimer("W high mode");
        }

        startTimer("V I/O");
        A2AVectorsIo::writeElement(par().output + "_v", vBin, b, vm().getTrajectory());
        stopTimer("V I/O");

        startTimer("W I/O");
        A2AVectorsIo::writeElement(par().output + "_w", wBin, b, vm().getTrajectory());
        stopTimer("W I/O");
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSolver_A2AVectorsCoarseHigh_hpp_
