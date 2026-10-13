/*
 * BubbleContractorDense.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#include <Hadrons/Global.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Hadrons/DiskVector.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/TimerArray.hpp>

using namespace Grid;
using namespace Hadrons;

#define TIME_MOD(t) (((t) + par.global.nt) % par.global.nt)

namespace BubbleContractorDense
{
    const std::size_t nSc = Ns*Nc;

    class TrajRange: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(TrajRange,
                                        unsigned int, start,
                                        unsigned int, end,
                                        unsigned int, step);
    };

    class GlobalPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(GlobalPar,
                                        TrajRange, trajCounter,
                                        unsigned int, nt,
                                        unsigned int, nHit,
                                        std::string, diskVectorDir,
                                        std::string, output);
    };

    class A2AMatrixPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMatrixPar,
                                        std::string, file,
                                        std::string, dataset,
                                        unsigned int, cacheSize,
                                        bool, timeSliceIO,
                                        unsigned int, nLow,
                                        std::string, name);
    };

    class ProductPar: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(ProductPar,
                                        std::string, terms,
                                        std::vector<std::string>, times);
    };

    class CorrelatorResult: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(CorrelatorResult,
                                        std::vector<BubbleContractorDense::A2AMatrixPar>,  a2aMatrix,
                                        ProductPar, contraction,
                                        std::vector<unsigned int>, times,
                                        std::vector<ComplexD>, correlator);
    };
}

struct BubbleContractorDensePar
{
    BubbleContractorDense::GlobalPar                  global;
    std::vector<BubbleContractorDense::A2AMatrixPar>  a2aMatrix;
    std::vector<BubbleContractorDense::ProductPar>    product;
};

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq,
                 const std::vector<std::set<unsigned int>> &times,
                 std::vector<unsigned int> &current,
                 const unsigned int depth)
{
    if (depth > 0)
    {
        for (auto t: times[times.size() - depth])
        {
            current[times.size() - depth] = t;
            makeTimeSeq(timeSeq, times, current, depth - 1);
        }
    }
    else
    {
        timeSeq.push_back(current);
    }
}

void makeTimeSeq(std::vector<std::vector<unsigned int>> &timeSeq,
                 const std::vector<std::set<unsigned int>> &times)
{
    std::vector<unsigned int> current(times.size());

    makeTimeSeq(timeSeq, times, current, times.size());
}

void saveCorrelator(const BubbleContractorDense::CorrelatorResult &result,
                    const std::string dir, const std::string fileStem,
                    const unsigned int traj)
{
    std::string filename;

    filename = dir + "/" + ModuleBase::resultFilename(fileStem, traj);
    std::cout << "Saving correlator to '" << filename << "'" << std::endl;
    makeFileDir(dir);
    ResultWriter writer(filename);
    write(writer, fileStem, result);
}

std::set<unsigned int> parseTimeRange(const std::string str, const unsigned int nt)
{
    std::regex               rex("([0-9]+)|(([0-9]+)\\.\\.([0-9]+))");
    std::smatch              sm;
    std::vector<std::string> rstr = strToVec<std::string>(str);
    std::set<unsigned int>   tSet;

    for (auto &s: rstr)
    {
        std::regex_match(s, sm, rex);
        if (sm[1].matched)
        {
            unsigned int t;

            t = std::stoi(sm[1].str());
            if (t >= nt)
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            tSet.insert(t);
        }
        else if (sm[2].matched)
        {
            unsigned int ta, tb;

            ta = std::stoi(sm[3].str());
            tb = std::stoi(sm[4].str());
            if ((ta >= nt) or (tb >= nt))
            {
                HADRONS_ERROR(Range, "time out of range (from expression '" + str + "')");
            }
            for (unsigned int ti = ta; ti <= tb; ++ti)
            {
                tSet.insert(ti);
            }
        }
    }

    return tSet;
}

std::string sliceFilename(const std::string stem, const unsigned int t)
{
    const std::string  ext = ".h5";
    std::ostringstream ss;

    if ((stem.size() < ext.size())
        or (stem.compare(stem.size() - ext.size(), ext.size(), ext) != 0))
    {
        HADRONS_ERROR(Io, "timeSliceIO file '" + stem + "' does not end in '.h5'");
    }
    ss << stem.substr(0, stem.size() - ext.size()) << ".t"
       << std::setfill('0') << std::setw(4) << t << ext;

    return ss.str();
}

struct Sec
{
    Sec(const double usec)
    {
        seconds = usec/1.0e6;
    }

    double seconds;
};

inline std::ostream & operator<< (std::ostream& s, const Sec &&sec)
{
    s << std::setw(10) << sec.seconds << " sec";

    return s;
}

struct Flops
{
    Flops(const double flops, const double fusec)
    {
        gFlopsPerSec = flops/fusec/1.0e3;
    }

    double gFlopsPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Flops &&f)
{
    s << std::setw(10) << f.gFlopsPerSec << " GFlop/s";

    return s;
}

struct Bytes
{
    Bytes(const double bytes, const double busec)
    {
        gBytesPerSec = bytes/busec*1.0e6/1024/1024/1024;
    }

    double gBytesPerSec;
};

inline std::ostream & operator<< (std::ostream& s, const Bytes &&b)
{
    s << std::setw(10) << b.gBytesPerSec << " GB/s";

    return s;
}

int main(int argc, char* argv[])
{
    std::string parFilename;

    if (argc != 2)
    {
        std::cerr << "usage: " << argv[0] << " <parameter file>";
        std::cerr << std::endl;

        return EXIT_FAILURE;
    }
    parFilename = argv[1];

    BubbleContractorDensePar par;
    XmlReader                reader(parFilename);

    read(reader, "global",    par.global);
    read(reader, "a2aMatrix", par.a2aMatrix);
    read(reader, "product",   par.product);

    std::map<std::string, EigenDiskVector<ComplexD>> a2aMat;
    std::map<std::string, A2AModeSpace>              modeSpace;

    for (auto &p: par.a2aMatrix)
    {
        std::string dirName = par.global.diskVectorDir + "/" + p.name;
        A2AModeSpace space;

        space.nLow = p.nLow;
        space.nSc  = BubbleContractorDense::nSc;
        space.nt   = par.global.nt;
        space.nHit = par.global.nHit;
        a2aMat.emplace(p.name, EigenDiskVector<ComplexD>(dirName, par.global.nt, p.cacheSize));
        modeSpace.emplace(p.name, space);
        std::cout << "A2A matrix '" << p.name << "': dense " << space.denseSize()
                  << ", expanded " << space.expandedSize() << std::endl;
    }

    for (unsigned int traj = par.global.trajCounter.start;
         traj < par.global.trajCounter.end; traj += par.global.trajCounter.step)
    {
        std::cout << ":::::::: Trajectory " << traj << std::endl;

        for (auto &p: par.a2aMatrix)
        {
            std::string filename = p.file;
            double      t, tTot = 0.;
            size_t      bytes = 0;

            tokenReplace(filename, "traj", traj);
            std::cout << "======== Loading '" << filename << "'" << std::endl;

            if (p.timeSliceIO)
            {
                EigenDiskVector<ComplexD>::Matrix buf;
                unsigned int                      ni = 0, nj = 0;

                for (unsigned int tp1 = par.global.nt; tp1 > 0; --tp1)
                {
                    unsigned int tt = tp1 - 1;

                    A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(sliceFilename(filename, tt),
                                                            p.dataset, 1, ni, nj);

                    a2aIo.loadSlice(buf, &t);
                    ni     = a2aIo.getNi();
                    nj     = a2aIo.getNj();
                    bytes += a2aIo.getSize();
                    tTot  += t;
                    a2aMat.at(p.name)[tt] = buf;
                }
            }
            else
            {
                A2AMatrixIo<HADRONS_A2AM_IO_TYPE> a2aIo(filename, p.dataset, par.global.nt);

                a2aIo.load(a2aMat.at(p.name), &t);
                bytes = a2aIo.getSize();
                tTot  = t;
            }
            std::cout << "Read " << bytes << " bytes in " << tTot/1.0e6
                    << " sec, " << bytes/tTot*1.0e6/1024/1024 << " MB/s" << std::endl;
        }

        for (auto &p: par.product)
        {
            std::vector<std::string>               term = strToVec<std::string>(p.terms);
            std::vector<std::set<unsigned int>>    times;
            std::vector<std::vector<unsigned int>> timeSeq;
            std::vector<A2AMatrixTr<ComplexD>>     lastTerm(par.global.nt);
            A2AMatrix<ComplexD>                    prod;
            TimerArray                             tAr;
            double                                 fusec, busec, flops, bytes;
            BubbleContractorDense::CorrelatorResult result;

            tAr.startTimer("Total");
            std::cout << "======== Bubble tr(";
            for (unsigned int g = 0; g < term.size(); ++g)
            {
                std::cout << term[g] << ((g == term.size() - 1) ? ')' : '*');
            }
            std::cout << std::endl;
            if ((term.size() < 1) or (term.size() > 2))
            {
                HADRONS_ERROR(Size, "BubbleContractorDense supports 1 or 2 terms, got "
                            + std::to_string(term.size())
                            + " (more needs a dense matrix product)");
            }
            if (term.size() != p.times.size() + 1)
            {
                HADRONS_ERROR(Size, "number of terms (" + std::to_string(term.size())
                            + ") different from number of times ("
                            + std::to_string(p.times.size() + 1) + ")");
            }
            for (auto &s: term)
            {
                if (modeSpace.find(s) == modeSpace.end())
                {
                    HADRONS_ERROR(Definition, "no a2aMatrix named '" + s + "'");
                }
            }

            std::string fileStem = "";
            for (unsigned int i = 0; i < term.size() - 1; i++)
            {
                fileStem += term[i] + "_" + p.times[i] + "_";
            }
            fileStem += term.back();

            const A2AModeSpace &spaceI = modeSpace.at(term[0]);
            const A2AModeSpace &spaceJ = modeSpace.at(term.back());

            for (auto &s: p.times)
            {
                times.push_back(parseTimeRange(s, par.global.nt));
            }
            for (auto &m: par.a2aMatrix)
            {
                if (std::find(result.a2aMatrix.begin(), result.a2aMatrix.end(), m) == result.a2aMatrix.end())
                {
                    result.a2aMatrix.push_back(m);
                    tokenReplace(result.a2aMatrix.back().file, "traj", traj);
                }
            }
            result.contraction = p;
            result.correlator.resize(par.global.nt, 0.);

            makeTimeSeq(timeSeq, times);

            if (term.size() > 1)
            {
                std::cout << "* Caching transposed last term" << std::endl;
                for (unsigned int t = 0; t < par.global.nt; ++t)
                {
                    tAr.startTimer("Disk vector overhead");
                    const A2AMatrix<ComplexD> &ref = a2aMat.at(term.back())[t];
                    tAr.stopTimer("Disk vector overhead");

                    tAr.startTimer("Transpose caching");
                    lastTerm[t].resize(ref.rows(), ref.cols());
                    thread_for( j,ref.cols(),{
                      for (unsigned int i = 0; i < ref.rows(); ++i)
                      {
                          lastTerm[t](i, j) = ref(i, j);
                      }
		    });
                    tAr.stopTimer("Transpose caching");
                }
                bytes = par.global.nt*lastTerm[0].rows()*lastTerm[0].cols()*sizeof(ComplexD);
                std::cout << Sec(tAr.getDTimer("Transpose caching")) << " "
                          << Bytes(bytes, tAr.getDTimer("Transpose caching")) << std::endl;
            }

            for (unsigned int i = 0; i < timeSeq.size(); ++i)
            {
                auto &t = timeSeq[i];

                result.times = t;
                for (unsigned int dt = 0; dt < par.global.nt; ++dt)
                {
                    result.correlator[dt] = 0.;
                }
                flops = 0.;
                bytes = 0.;
                fusec = tAr.getDTimer("tr");
                busec = tAr.getDTimer("tr");
                tAr.startTimer("Linear algebra");
                for (unsigned int dt = 0; dt < par.global.nt; ++dt)
                {
                    tAr.startTimer("Disk vector overhead");
                    prod = a2aMat.at(term[0])[TIME_MOD(dt)];
                    tAr.stopTimer("Disk vector overhead");

                    tAr.startTimer("tr");
                    if (term.size() == 1)
                    {
                        A2AContraction::accTrDense(result.correlator[dt], prod,
                                                   TIME_MOD(dt), spaceI);
                    }
                    else
                    {
                        unsigned int tLast = TIME_MOD(t[0] + dt);

                        A2AContraction::accTrMulDense(result.correlator[dt],
                                                      prod, TIME_MOD(dt),
                                                      lastTerm[tLast], tLast,
                                                      spaceI, spaceJ);
                        flops += A2AContraction::accTrMulDenseFlops(spaceI, spaceJ);
                        bytes += 2.*prod.rows()*prod.cols()*sizeof(ComplexD);
                    }
                    tAr.stopTimer("tr");
                }
                tAr.stopTimer("Linear algebra");
                if (term.size() > 1)
                {
                    std::cout << std::setw(8) << "traces"
                            << Sec(tAr.getDTimer("tr") - busec) << " "
                            << Flops(flops, tAr.getDTimer("tr") - fusec) << " "
                            << Bytes(bytes, tAr.getDTimer("tr") - busec) << std::endl;
                }
                saveCorrelator(result, par.global.output, fileStem, traj);
            }
            tAr.stopTimer("Total");
            printTimeProfile(tAr.getTimings(), tAr.getTimer("Total"));
        }
    }

    return EXIT_SUCCESS;
}
