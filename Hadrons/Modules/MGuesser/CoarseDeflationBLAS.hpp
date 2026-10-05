/*
 * CoarseDeflationBLAS.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
 * Author: Raoul Hodgson <raoul.hodgson@ed.ac.uk>
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
#ifndef Hadrons_MGuesser_CoarseDeflationBLAS_hpp_
#define Hadrons_MGuesser_CoarseDeflationBLAS_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Coarse deflation guesser, BLAS coefficients                *
 ******************************************************************************/
/*
 * Same guess as MGuesser::CoarseDeflation, with the coarse-space projection
 * done by GEMM instead of a loop of inner products. MGuesser::CoarseDeflation
 * accumulates sum_i evec_i (evec_i, src)/eval_i as nev separate innerProduct
 * calls per source, and every innerProduct ends in its own global reduction, so
 * a solve costs nev reductions whether or not the sources are batched.
 * Grid's MultiRHSDeflation instead forms the whole coefficient matrix
 * C = E^dag R with one GEMM, reduces it with a single GlobalSumVector over all
 * nev x nrhs entries, scales by 1/eval on the host, and reconstructs
 * G = E C with a second GEMM. One reduction per call replaces nev per source,
 * which is where the time goes (measured: axpy is ~80% of the guesser).
 *
 * Block projection and promotion are left as they are, on the fine subspace,
 * because they are a small part of the cost. Moving them to
 * Grid's MultiRHSBlockProject is the remaining step and belongs in this module:
 * it would also let the Lattice eigenpack be released once the BLAS buffers are
 * packed, since nothing here would reference it any more.
 *
 * Cost of the coefficient path is one extra copy of the coarse eigenbasis, as
 * pinned device memory (BLAS_E, nev x coarse_vol x nbasis complex). With a
 * blockSize that does not block Ls the coarse grid keeps the full fifth
 * dimension, so this is not a small buffer -- size it before running.
 *
 * Numerically this is NOT bit-compatible with MGuesser::CoarseDeflation: a GEMM
 * and a loop of axpy sum in different orders. Compare final residuals and
 * iteration counts, not checksums.
 */
BEGIN_MODULE_NAMESPACE(MGuesser)

class CoarseDeflationBLASPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(CoarseDeflationBLASPar,
                                    std::string, eigenPack,
                                    unsigned int, size);
};

template <typename EPack>
class TCoarseDeflationBLAS: public Module<CoarseDeflationBLASPar>
{
public:
    typedef typename EPack::Field Field;
    typedef typename EPack::CoarseField CoarseField;
public:
    // constructor
    TCoarseDeflationBLAS(const std::string name);
    // destructor
    virtual ~TCoarseDeflationBLAS(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual DependencyMap getObjectDependencies(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(CoarseDeflationBLAS200 , ARG(TCoarseDeflationBLAS<CoarseFermionEigenPack<FIMPL ,200>>), MGuesser);
MODULE_REGISTER_TMP(CoarseDeflationBLAS200F, ARG(TCoarseDeflationBLAS<CoarseFermionEigenPack<FIMPLF,200>>), MGuesser);

/******************************************************************************
 *                            The guesser itself                              *
 ******************************************************************************/
template<class FineField, class CoarseField>
class LocalCoherenceDeflatedGuesserBLAS: public LinearFunction<FineField> {
private:
    const std::vector<FineField>   &subspace_;
    const unsigned int             epackSize_;
    GridBase                       *coarseGrid_;
    MultiRHSDeflation<CoarseField> deflate_;

public:
    using LinearFunction<FineField>::operator();

    LocalCoherenceDeflatedGuesserBLAS(const std::vector<FineField> &subspace,
                                      const std::vector<CoarseField> &evec_coarse,
                                      const std::vector<RealD> &eval)
    : LocalCoherenceDeflatedGuesserBLAS(subspace,evec_coarse,eval,evec_coarse.size())
    {}

    //////////////////////////////////////////////////////////////////////////
    // The coarse eigenbasis is packed into BLAS_E here, once. The eigenpack
    // has been read by the time a guesser is constructed: the virtual machine
    // runs one module at a time, setup then execute, so the pack's loader has
    // already executed. Only subspace_ is kept by reference afterwards, for the
    // block projection; no reference to the coarse eigenvectors survives.
    //
    // Allocate + ImportEigenVector rather than ImportEigenBasis: the latter
    // logs one line per eigenvector through a bare std::cout, with no rank
    // filtering, which at production mode counts is hundreds of thousands of
    // lines.
    //////////////////////////////////////////////////////////////////////////
    LocalCoherenceDeflatedGuesserBLAS(const std::vector<FineField> &subspace,
                                      const std::vector<CoarseField> &evec_coarse,
                                      const std::vector<RealD> &eval,
                                      unsigned int epackSize)
    : subspace_(subspace), epackSize_(epackSize)
    {
        assert(evec_coarse.size() == eval.size());
        assert(epackSize_ <= evec_coarse.size());

        coarseGrid_ = evec_coarse[0].Grid();

        auto &evec = const_cast<std::vector<CoarseField>&>(evec_coarse);
        auto &ev   = const_cast<std::vector<RealD>&>(eval);
        double t   = -usecond();

        deflate_.Allocate(epackSize_, coarseGrid_);
        for (unsigned int i = 0; i < epackSize_; ++i)
        {
            deflate_.ImportEigenVector(evec[i], ev[i], i);
        }
        t += usecond();
        LOG(Message) << "LocalCoherenceDeflatedGuesserBLAS: packed " << epackSize_
                     << " coarse eigenvectors in " << t/1.0e6 << " s" << std::endl;
    }

    virtual void operator()(const FineField &src,FineField &guess) {
        std::vector<FineField> srcVec   = {src};
        std::vector<FineField> guessVec = {guess};

        (*this)(srcVec,guessVec);

        guess = guessVec[0];
    }

    virtual void operator() (const std::vector<FineField> &src, std::vector<FineField> &guess)
    {
        assert(src.size() == guess.size());

        unsigned int sourceSize = src.size();

        double time_project = 0.;
        double time_coeff   = 0.;
        double time_promote = 0.;

        std::vector<CoarseField> src_coarse;   src_coarse.reserve(sourceSize);
        std::vector<CoarseField> guess_coarse; guess_coarse.reserve(sourceSize);
        for (unsigned int k = 0; k < sourceSize; ++k) {
            src_coarse.emplace_back(coarseGrid_);
            guess_coarse.emplace_back(coarseGrid_);
        }

        time_project -= usecond();
        batchBlockProject(src_coarse,src,subspace_);
        time_project += usecond();

        // DeflateSources overwrites every site of each guess, so guess_coarse
        // needs no initialisation.
        time_coeff -= usecond();
        deflate_.DeflateSources(src_coarse,guess_coarse);
        time_coeff += usecond();

        time_promote -= usecond();
        batchBlockPromote(guess_coarse,guess,subspace_);
        time_promote += usecond();

        // The BLAS path moves raw device memory and never sets a checkerboard
        for (unsigned int k = 0; k < sourceSize; ++k)
            guess[k].Checkerboard() = src[k].Checkerboard();

        LOG(Message) << "LocalCoherenceDeflatedGuesserBLAS: Total projection   time " << time_project/1.e6 << " s" <<  std::endl;
        LOG(Message) << "LocalCoherenceDeflatedGuesserBLAS: Total coefficient  time " << time_coeff/1.e6   << " s" <<  std::endl;
        LOG(Message) << "LocalCoherenceDeflatedGuesserBLAS: Total promotion    time " << time_promote/1.e6 << " s" <<  std::endl;
    }
};

/******************************************************************************
 *                 TCoarseDeflationBLAS implementation                        *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename EPack>
TCoarseDeflationBLAS<EPack>::TCoarseDeflationBLAS(const std::string name)
: Module<CoarseDeflationBLASPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename EPack>
std::vector<std::string> TCoarseDeflationBLAS<EPack>::getInput(void)
{
    std::vector<std::string> in = {par().eigenPack};

    return in;
}

template <typename EPack>
std::vector<std::string> TCoarseDeflationBLAS<EPack>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename EPack>
DependencyMap TCoarseDeflationBLAS<EPack>::getObjectDependencies(void)
{
    DependencyMap dep;

    dep.insert({par().eigenPack, getName()});

    return dep;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename EPack>
void TCoarseDeflationBLAS<EPack>::setup(void)
{
    LOG(Message) << "Setting up local coherence deflation guesser (BLAS coefficients) "
                 << "with eigenpack '" << par().eigenPack << "'" << std::endl;

    auto &epack = envGet(EPack, par().eigenPack);
    envCreateDerived(LinearFunction<Field>, ARG(LocalCoherenceDeflatedGuesserBLAS<Field, CoarseField>), getName(),
                     env().getObjectLs(par().eigenPack), epack.evec, epack.evecCoarse, epack.evalCoarse);

}

// execution ///////////////////////////////////////////////////////////////////
template <typename EPack>
void TCoarseDeflationBLAS<EPack>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MGuesser_CoarseDeflationBLAS_hpp_
