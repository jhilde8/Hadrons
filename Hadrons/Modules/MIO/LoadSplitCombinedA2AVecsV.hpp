/*
 * LoadSplitCombinedA2AVecsV.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2026
 *
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
#ifndef Hadrons_MIO_LoadSplitCombinedA2AVecsV_hpp_
#define Hadrons_MIO_LoadSplitCombinedA2AVecsV_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AVectors.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  MIO::LoadCombinedA2AVecsV with the reads done on split-grid partitions.  *
 *  Same parameters, files, layout and 1/nHit normalization, plus ioSplit,   *
 *  the partition processor layout (e.g. "4 4 4 4").                         *
 *                                                                            *
 *  A per-file read is bound by a fixed cost rather than by the reader's     *
 *  size, so reading many files at once raises the aggregate rate. Each      *
 *  partition reads a different file into one split-grid field, one file per *
 *  partition per batch, and Grid_unsplit moves the batch into its full-grid *
 *  slots -- each byte read crosses the network once more, in whole          *
 *  full-grid local blocks. Nothing is split: the data starts on the split   *
 *  grid, so the unsplit is the only redistribution.                         *
 *                                                                            *
 *  Timers: "read" is this rank's partition's reads; "read wait" is a        *
 *  barrier before each unsplit, so it holds the slowest partition's read    *
 *  less this one; "unsplit" is then the redistribution alone.               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MIO)

class LoadSplitCombinedA2AVecsVPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(LoadSplitCombinedA2AVecsVPar,
                                    std::string,               lowFilestem,
                                    unsigned int,              nLow,
                                    std::string,               highStem,
                                    std::vector<std::string>, highExtensions,
                                    unsigned int,              highSize,
                                    unsigned int,              nHit,
                                    std::string,               ioSplit);
};

template <typename FImpl>
class TLoadSplitCombinedA2AVecsV: public Module<LoadSplitCombinedA2AVecsVPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TLoadSplitCombinedA2AVecsV(const std::string name);
    // destructor
    virtual ~TLoadSplitCombinedA2AVecsV(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    void readBlock(const std::string &stem, std::vector<FermionField> &out,
                   const unsigned int offset, const unsigned int size);
};

MODULE_REGISTER_TMP(LoadSplitCombinedA2AVecsV, ARG(TLoadSplitCombinedA2AVecsV<FIMPL>), MIO);

/******************************************************************************
 *                TLoadSplitCombinedA2AVecsV implementation                   *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TLoadSplitCombinedA2AVecsV<FImpl>::TLoadSplitCombinedA2AVecsV(const std::string name)
: Module<LoadSplitCombinedA2AVecsVPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TLoadSplitCombinedA2AVecsV<FImpl>::getInput(void)
{
    std::vector<std::string> in;

    return in;
}

template <typename FImpl>
std::vector<std::string> TLoadSplitCombinedA2AVecsV<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TLoadSplitCombinedA2AVecsV<FImpl>::setup(void)
{
    unsigned int total = par().nLow + par().highExtensions.size() * par().highSize;
    auto         grid  = envGetGrid(FermionField);

    envCreate(std::vector<FermionField>, getName(), 1, 0, grid);
    auto &out = envGet(std::vector<FermionField>, getName());
    out.reserve(total);
    for (unsigned int i = 0; i < total; ++i)
    {
        out.emplace_back(grid, CpuWrite);
    }
}

// execution ///////////////////////////////////////////////////////////////////
// Files 0..size-1 of stem into out[offset..offset+size). Batch b reads file
// b + part on partition part; partitions past the end idle but still join the
// collective unsplit. swap hands the batch slot the out slot's unused
// allocation, so nothing is allocated per batch.
template <typename FImpl>
void TLoadSplitCombinedA2AVecsV<FImpl>::readBlock(const std::string &stem,
                                                  std::vector<FermionField> &out,
                                                  const unsigned int offset,
                                                  const unsigned int size)
{
    int           traj = vm().getTrajectory();
    auto          grid = envGetGrid(FermionField);
    GridCartesian splitGrid(grid->FullDimensions(), grid->_simd_layout,
                            strToVec<int>(par().ioSplit), *grid);
    unsigned int  nPart = grid->ProcessorCount()/splitGrid.ProcessorCount();
    unsigned int  part  = GridSplitVectorIndex(grid, &splitGrid);
    FermionField  split(&splitGrid);
    std::vector<FermionField> batch(nPart, grid);

    LOG(Message) << "Split read: " << nPart << " partitions of layout "
                 << par().ioSplit << std::endl;
    for (unsigned int b = 0; b < size; b += nPart)
    {
        startTimer("read");
        if (b + part < size)
        {
            A2AVectorsIo::readElement(stem, split, b + part, traj);
        }
        stopTimer("read");
        startTimer("read wait");
        grid->Barrier();
        stopTimer("read wait");
        startTimer("unsplit");
        Grid_unsplit(batch, split);
        stopTimer("unsplit");
        for (unsigned int p = 0; p < nPart && b + p < size; ++p)
        {
            std::swap(out[offset + b + p], batch[p]);
        }
    }
}

template <typename FImpl>
void TLoadSplitCombinedA2AVecsV<FImpl>::execute(void)
{
    auto &out = envGet(std::vector<FermionField>, getName());
    unsigned int offset = 0;

    if (par().nLow > 0)
    {
        LOG(Message) << "Loading " << par().nLow << " low-mode A2A vectors from '"
                     << par().lowFilestem << "'" << std::endl;
        readBlock(par().lowFilestem, out, offset, par().nLow);
        offset += par().nLow;
    }

    Real norm = (par().nHit > 0) ? (1.0/par().nHit) : 1.0;

    for (unsigned int h = 0; h < par().highExtensions.size(); ++h)
    {
        std::string filestem = par().highStem + par().highExtensions[h];

        LOG(Message) << "Loading " << par().highSize << " high-mode A2A vectors from '"
                     << filestem << "'" << std::endl;
        readBlock(filestem, out, offset, par().highSize);
        if (norm != 1.0)
        {
            LOG(Message) << "Applying 1/nHit = " << norm
                         << " hit-average normalization" << std::endl;
            for (unsigned int i = offset; i < offset + par().highSize; ++i)
            {
                out[i] = norm*out[i];
            }
        }
        offset += par().highSize;
    }
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MIO_LoadSplitCombinedA2AVecsV_hpp_
