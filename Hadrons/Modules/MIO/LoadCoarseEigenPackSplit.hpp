/*
 * LoadCoarseEigenPackSplit.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MIO_LoadCoarseEigenPackSplit_hpp_
#define Hadrons_MIO_LoadCoarseEigenPackSplit_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/EigenPack.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  MIO::LoadCoarseEigenPack with the coarse vectors read on split-grid        *
 *  partitions.                                                                *
 *                                                                             *
 *  The coarse read is latency-bound: on the full grid every file (~210 MB on  *
 *  64I) is a collective read of a few hundred kB per rank, 2000 times in      *
 *  series. Here each partition of layout ioSplit (4D, e.g. "1 1 2 4" for 8    *
 *  ranks) reads whole files, nPart at a time, and Grid_unsplit hands them to   *
 *  the full grid; per-partition timings are logged as for the A2A vectors.    *
 *  The fine vectors are read as before.                                       *
 *                                                                             *
 *  Multi-file packs only, with I/O type equal to the vector type (the 200 and *
 *  200F registrations). Kept apart from MIO::LoadCoarseEigenPack until        *
 *  tested; Test_eigenpack_split_read checks it bitwise against readCoarse.    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadCoarseEigenPackSplitPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadCoarseEigenPackSplitPar,
                                    std::string,  filestem,
                                    unsigned int, sizeFine,
                                    unsigned int, sizeCoarse,
                                    bool,         redBlack,
                                    unsigned int, Ls,
                                    std::string,  blockSize,
                                    bool,         orthogonalise,
                                    std::string,  gaugeXform,
                                    std::string,  ioSplit);
};

template <typename Pack, typename GImpl>
class TLoadCoarseEigenPackSplit: public Module<LoadCoarseEigenPackSplitPar>
{
public:
    typedef typename Pack::Field                Field;
    typedef typename Pack::CoarseField          CoarseField;
    typedef CoarseEigenPack<Field, CoarseField> BasePack;
    template <typename vtype>
    using iImplScalar = iScalar<iScalar<iScalar<vtype>>>;
    typedef iImplScalar<typename Pack::Field::vector_type> SiteComplex;

    GAUGE_TYPE_ALIASES(GImpl, );
    typedef typename GImpl::GaugeLinkField GaugeMat;
public:
    // constructor
    TLoadCoarseEigenPackSplit(const std::string name);
    // destructor
    virtual ~TLoadCoarseEigenPackSplit(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackSplit200,  ARG(TLoadCoarseEigenPackSplit<CoarseFermionEigenPack<FIMPL, 200>, GIMPL>), MIO);
MODULE_REGISTER_TMP(LoadCoarseFermionEigenPackSplit200F, ARG(TLoadCoarseEigenPackSplit<CoarseFermionEigenPack<FIMPLF, 200>, GIMPLF>), MIO);

/******************************************************************************
 *                 TLoadCoarseEigenPackSplit implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
TLoadCoarseEigenPackSplit<Pack,GImpl>::TLoadCoarseEigenPackSplit(const std::string name)
: Module<LoadCoarseEigenPackSplitPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
std::vector<std::string> TLoadCoarseEigenPackSplit<Pack,GImpl>::getInput(void)
{
    std::vector<std::string> in;

    if (!par().gaugeXform.empty())
    {
        in = {par().gaugeXform};
    }

    return in;
}

template <typename Pack, typename GImpl>
std::vector<std::string> TLoadCoarseEigenPackSplit<Pack,GImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadCoarseEigenPackSplit<Pack,GImpl>::setup(void)
{
    auto blockSize = strToVec<int>(par().blockSize);

    if (par().ioSplit.empty())
    {
        HADRONS_ERROR(Argument, "ioSplit is empty, use MIO::LoadCoarseEigenPack");
    }
    envCreateDerived(BasePack, Pack, getName(), par().Ls, par().sizeFine,
                     par().sizeCoarse, envGetRbGrid(Field, par().Ls),
                     envGetCoarseGrid(CoarseField, blockSize, par().Ls));

    GridBase* grid   = envGetGrid(Field, par().Ls);
    GridBase* gridRb = envGetRbGrid(Field, par().Ls);
    if (!par().gaugeXform.empty())
    {
        envTmp(GaugeMat,    "tmpXform", par().Ls, grid);
        if (par().redBlack)
        {
            envTmp(GaugeMat, "tmpXformOdd", par().Ls, gridRb);
        }
    }
}

// execution ///////////////////////////////////////////////////////////////////
template <typename Pack, typename GImpl>
void TLoadCoarseEigenPackSplit<Pack,GImpl>::execute(void)
{
    auto                 blockSize = strToVec<int>(par().blockSize);
    auto                 cg        = envGetCoarseGrid(CoarseField, blockSize, par().Ls);
    auto                 cg4d      = envGetCoarseGrid(CoarseField, blockSize);
    auto                 &epack    = envGetDerived(BasePack, Pack, getName());
    Lattice<SiteComplex> dummy(cg);
    std::vector<double>  t;

    GridCartesian                  split4d(cg4d->FullDimensions(), cg4d->_simd_layout,
                                           strToVec<int>(par().ioSplit), *cg4d);
    std::unique_ptr<GridCartesian> split(SpaceTimeGrid::makeFiveDimGrid(
                                             cg->FullDimensions()[0], &split4d));
    unsigned int                   nPart = cg->ProcessorCount()/split->ProcessorCount();
    double                         bytes = static_cast<double>(cg->gSites())
                                           *sizeof(typename CoarseField::scalar_object);

    LOG(Message) << "Split coarse read: " << nPart << " partitions of layout "
                 << par().ioSplit << std::endl;

    startTimer("Read fine");
    epack.readFine(par().filestem, true, vm().getTrajectory());
    stopTimer("Read fine");
    startTimer("Read coarse");
    epack.readCoarseSplit(par().filestem, split.get(), t, vm().getTrajectory(), this);
    stopTimer("Read coarse");
    A2AVectorsIo::logPartitions("Read coarse", 0, par().sizeCoarse, nPart, t, bytes);

    if (!par().gaugeXform.empty())
    {
        LOG(Message) << "Applying gauge transformation to fine eigenvectors " << getName()
                     << " using " << par().gaugeXform << std::endl;

        auto &xform = envGet(GaugeMat, par().gaugeXform);
        epack.gaugeTransform(xform);
    }

    if (par().orthogonalise) {
        LOG(Message) << "Block Gramm-Schmidt pass 1"<< std::endl;
        blockOrthogonalise(dummy, epack.evec);
        LOG(Message) << "Block Gramm-Schmidt pass 2"<< std::endl;
        blockOrthogonalise(dummy, epack.evec);
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadCoarseEigenPackSplit_hpp_
