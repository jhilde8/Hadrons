/*
 * A2ANewMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: ferben <ferben@debian.felix.com>
 * Author: paboyle <paboyle@ph.ed.ac.uk>
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
#ifndef Hadrons_MContraction_A2ANewMesonField_hpp_
#define Hadrons_MContraction_A2ANewMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>
#include <Grid/algorithms/blas/A2ASpatialSum.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                     All-to-all meson field creation (new, no A2AMatrixBlockComputation)
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2ANewMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewMesonFieldPar,
                                    int,                     block,
                                    std::string,             left,
                                    std::string,             right,
                                    std::string,             output,
                                    std::string,             gammas,
                                    std::vector<std::string>, mom);
};

class A2ANewMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2ANewMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma);
};

template <typename FImpl>
class TA2ANewMesonField : public Module<A2ANewMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TA2ANewMesonField(const std::string name);
    virtual ~TA2ANewMesonField(void){};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
private:
    bool                                      hasPhase_{false};
    std::string                               momphName_;
    std::vector<Gamma::Algebra>               gamma_;
    std::vector<std::vector<Real>>            mom_;
    A2ASpatialSum<typename FImpl::SiteSpinor> spatial_sum_;
};

MODULE_REGISTER(A2ANewMesonField, ARG(TA2ANewMesonField<FIMPL>), MContraction);

/******************************************************************************
 *                  TA2ANewMesonField implementation                          *
 ******************************************************************************/
template <typename FImpl>
TA2ANewMesonField<FImpl>::TA2ANewMesonField(const std::string name)
: Module<A2ANewMesonFieldPar>(name)
, momphName_(name + "_momph")
{}

template <typename FImpl>
std::vector<std::string> TA2ANewMesonField<FImpl>::getInput(void)
{
    return {par().left, par().right};
}

template <typename FImpl>
std::vector<std::string> TA2ANewMesonField<FImpl>::getOutput(void)
{
    return {};
}

template <typename FImpl>
void TA2ANewMesonField<FImpl>::setup(void)
{
    gamma_.clear();
    mom_.clear();
    if (par().gammas == "all")
    {
        gamma_ = {
            Gamma::Algebra::Gamma5,
            Gamma::Algebra::Identity,
            Gamma::Algebra::GammaX,
            Gamma::Algebra::GammaY,
            Gamma::Algebra::GammaZ,
            Gamma::Algebra::GammaT,
            Gamma::Algebra::GammaXGamma5,
            Gamma::Algebra::GammaYGamma5,
            Gamma::Algebra::GammaZGamma5,
            Gamma::Algebra::GammaTGamma5,
            Gamma::Algebra::SigmaXY,
            Gamma::Algebra::SigmaXZ,
            Gamma::Algebra::SigmaXT,
            Gamma::Algebra::SigmaYZ,
            Gamma::Algebra::SigmaYT,
            Gamma::Algebra::SigmaZT
        };
    }
    else
    {
        gamma_ = strToVec<Gamma::Algebra>(par().gammas);
    }
    for (auto &pstr: par().mom)
    {
        auto p = strToVec<Real>(pstr);
        if (p.size() != env().getNd() - 1)
        {
            HADRONS_ERROR(Size, "Momentum has " + std::to_string(p.size())
                                + " components instead of "
                                + std::to_string(env().getNd() - 1));
        }
        mom_.push_back(p);
    }
    envCache(std::vector<ComplexField>, momphName_, 1,
             par().mom.size(), envGetGrid(ComplexField));
    envTmpLat(ComplexField, "coor");
}

template <typename FImpl>
void TA2ANewMesonField<FImpl>::execute(void)
{
    typedef typename FImpl::SiteSpinor      vobj;
    typedef typename vobj::vector_type      vector_type;
    typedef iSpinColourVector<vector_type>  SpinColourVector_v;
    typedef typename SpinColourVector_v::scalar_type scalar_t;

    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = envGetGrid(FermionField);

    int nt     = env().getDim().back();
    int N_i    = left.size();
    int N_j    = right.size();
    int ngamma = gamma_.size();
    int nmom   = mom_.size();
    int block  = par().block;

    LOG(Message) << "Computing all-to-all meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
        LOG(Message) << "  " << p << std::endl;
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
        LOG(Message) << "  " << g << std::endl;
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE))
                 << "/momentum/bilinear)" << std::endl;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (int j = 0; j < nmom; ++j)
        {
            Complex i(0.0, 1.0);
            envGetTmp(ComplexField, coor);
            ph[j] = Zero();
            for (unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        ComplexField last_abs_ph = ph[nmom - 1];
        for (int j = nmom - 1; j >= 1; --j)
            ph[j] = ph[j] * adj(ph[j - 1]); // apply the phase difference
        ph.push_back(adj(last_abs_ph)); // restore transition: undo final accumulated phase
        hasPhase_ = true;
        stopTimer("Momentum phases");
    }

    auto ionameFn = [this](const int m, const int g)
    {
        std::stringstream ss;
        ss << gamma_[g] << "_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const int m, const int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory())
               + "/" + ionameFn(m, g) + ".h5";
    };

    auto metadataFn = [this](const int m, const int g)
    {
        A2ANewMesonFieldMetadata md;
        for (auto pmu: mom_[m])
            md.momentum.push_back(pmu);
        md.gamma = gamma_[g];
        return md;
    };

    // Output buffer: one (nt, Nii, Njj) block at a time.
    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

    // Scratch right vectors for GammaRight output (zero-momentum base pack).
    std::vector<FermionField> gammaRight(block, grid);

    // Create output directory.
    std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());
    {
        std::string dummy = dirBase + "/mkdir.h5";
        makeFileDir(dummy, grid);
    }

    // Initialise one HDF5 file per (mom, gamma) pair before the block loops.
    for (int m = 0; m < nmom; m++)
    for (int g = 0; g < ngamma; g++)
    {
        std::string ioname   = ionameFn(m, g);
        std::string filename = filenameFn(m, g);
        A2ANewMesonFieldMetadata md = metadataFn(m, g);
#ifdef HADRONS_A2AM_PARALLEL_IO
        grid->Barrier();
        if (grid->ThisRank() == 0) {
#endif
        {
            A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j);
            io.initFile(md, block);
        }
#ifdef HADRONS_A2AM_PARALLEL_IO
        }
        grid->Barrier();
#endif
    }

    // Pre-pack flat phase arrays (one per momentum) for in-place buffer multiplication.
    // PackPhase mirrors PackVectors: SIMD/SIMT extraction -> one scalar per spatial site.
    startTimer("Pack phases");
    std::vector<deviceVector<scalar_t>> ph_flat(ph.size());
    for (int m = 0; m < (int)ph.size(); m++)
        A2ASpatialSum<SpinColourVector_v>::PackPhase(grid, ph[m], ph_flat[m]);
    stopTimer("Pack phases");

    // One-time allocation for the full block size; subsequent pointer rewrites are cheap.
    startTimer("Allocate");
    spatial_sum_.AllocateRight(block, grid);
    spatial_sum_.AllocateLeft(block);
    stopTimer("Allocate");

    // Loop order (jb, g, ib, m):
    //   AllocateRight + PackRight - once per (jb, g)
    //   AllocateLeft  + PackLeft  - once per ib
    //   Apply+GEMM+Restore        - once per (jb, g, ib, m)

    for (int jb = 0; jb < N_j; jb += block)
    {
        int Njj = std::min(N_j - jb, block);

        startTimer("Allocate");
        spatial_sum_.AllocateRight(Njj, grid);
        stopTimer("Allocate");

        for (int g = 0; g < ngamma; g++)
        {
            startTimer("GammaRight");
            for (int jj = 0; jj < Njj; jj++)
                A2Autils<FImpl>::GammaRight(gammaRight[jj], gamma_[g], right[jb + jj]);
            stopTimer("GammaRight");

            startTimer("Pack vectors");
            spatial_sum_.PackRight(gammaRight, 0, Njj);
            stopTimer("Pack vectors");

            for (int ib = 0; ib < N_i; ib += block)
            {
                int Nii = std::min(N_i - ib, block);

                startTimer("Allocate");
                spatial_sum_.AllocateLeft(Nii);
                stopTimer("Allocate");

                startTimer("Pack vectors");
                spatial_sum_.PackLeftConj(left, ib, Nii);
                stopTimer("Pack vectors");

                Eigen::Tensor<ComplexD, 3> block_result(nt, Nii, Njj);

                // ph_flat[0] is the absolute phase for the first momentum.
                // ph_flat[nmom] is the restore transition that steps LR_buf back to the
                // unphased state after the last Sum, so each ib block begins with clean
                // right vectors without needing a new PackRight call.
                startTimer("Phase");
                spatial_sum_.ApplyPhaseRight(ph_flat[0]);
                stopTimer("Phase");

                for (int m = 0; m < nmom; m++)
                {
                    startTimer("Sum");
                    spatial_sum_.Sum(block_result);
                    stopTimer("Sum");

                    startTimer("IO");
                    A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);
                    for (int t  = 0; t  < nt;  t++)
                    for (int ii = 0; ii < Nii; ii++)
                    for (int jj = 0; jj < Njj; jj++)
                        mf(0, 0, t, ii, jj) = block_result(t, ii, jj);

                    std::string ioname   = ionameFn(m, g);
                    std::string filename = filenameFn(m, g);

                    LOG(Message) << "MF block i=" << ib << " j=" << jb
                                 << " g=" << gamma_[g] << " m=" << m << std::endl;

#ifdef HADRONS_A2AM_PARALLEL_IO
                    grid->Barrier();
                    if (grid->ThisRank() == 0) {
#endif
                    {
                        A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j);
                        io.saveBlock(mf, 0, 0, ib, jb);
                    }
#ifdef HADRONS_A2AM_PARALLEL_IO
                    }
                    grid->Barrier();
#endif
                    stopTimer("IO");

                    startTimer("Phase");
                    spatial_sum_.ApplyPhaseRight(ph_flat[m + 1]);
                    stopTimer("Phase");
                } // m
            } // ib
        } // g
    } // jb
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

// =============================================================================
// OLD execute() — loop order (jb, g, m, ib) using PhaseContractRight.
// Kept for regression comparison.  Restore by swapping with the active execute()
// above and recompiling.
// =============================================================================
#if 0
template <typename FImpl>
void TA2ANewMesonField<FImpl>::execute_old(void)
{
    typedef typename FImpl::SiteSpinor      vobj;
    typedef typename vobj::vector_type      vector_type;
    typedef iSpinColourVector<vector_type>  SpinColourVector_v;

    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = envGetGrid(FermionField);

    int nt     = env().getDim().back();
    int N_i    = left.size();
    int N_j    = right.size();
    int ngamma = gamma_.size();
    int nmom   = mom_.size();
    int block  = par().block;

    auto &ph = envGet(std::vector<ComplexField>, momphName_);

    if (!hasPhase_)
    {
        startTimer("Momentum phases");
        for (int j = 0; j < nmom; ++j)
        {
            Complex i(0.0, 1.0);
            envGetTmp(ComplexField, coor);
            ph[j] = Zero();
            for (unsigned int mu = 0; mu < mom_[j].size(); mu++)
            {
                LatticeCoordinate(coor, mu);
                ph[j] = ph[j] + (mom_[j][mu]/env().getDim(mu))*coor;
            }
            ph[j] = exp((Real)(2*M_PI)*i*ph[j]);
        }
        hasPhase_ = true;
        stopTimer("Momentum phases");
    }

    auto ionameFn = [this](const int m, const int g)
    {
        std::stringstream ss;
        ss << gamma_[g] << "_";
        for (unsigned int mu = 0; mu < mom_[m].size(); ++mu)
            ss << mom_[m][mu] << ((mu == mom_[m].size() - 1) ? "" : "_");
        return ss.str();
    };

    auto filenameFn = [this, &ionameFn](const int m, const int g)
    {
        return par().output + "." + std::to_string(vm().getTrajectory())
               + "/" + ionameFn(m, g) + ".h5";
    };

    auto metadataFn = [this](const int m, const int g)
    {
        A2ANewMesonFieldMetadata md;
        for (auto pmu: mom_[m])
            md.momentum.push_back(pmu);
        md.gamma = gamma_[g];
        return md;
    };

    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

    std::vector<FermionField> phaseGammaRight(block, grid);

    std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());
    { std::string dummy = dirBase + "/mkdir.h5"; makeFileDir(dummy, grid); }

    for (int m = 0; m < nmom; m++)
    for (int g = 0; g < ngamma; g++)
    {
        std::string ioname   = ionameFn(m, g);
        std::string filename = filenameFn(m, g);
        A2ANewMesonFieldMetadata md = metadataFn(m, g);
#ifdef HADRONS_A2AM_PARALLEL_IO
        grid->Barrier();
        if (grid->ThisRank() == 0) {
#endif
        { A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j); io.initFile(md, block); }
#ifdef HADRONS_A2AM_PARALLEL_IO
        } grid->Barrier();
#endif
    }

    for (int jb = 0; jb < N_j; jb += block)
    {
        int Njj = std::min(N_j - jb, block);
        for (int g = 0; g < ngamma; g++)
        for (int m = 0; m < nmom; m++)
        {
            for (int jj = 0; jj < Njj; jj++)
                A2Autils<FImpl>::PhaseContractRight(
                    phaseGammaRight[jj], ph[m], gamma_[g], right[jb + jj]);
            std::string ioname   = ionameFn(m, g);
            std::string filename = filenameFn(m, g);
            for (int ib = 0; ib < N_i; ib += block)
            {
                int Nii = std::min(N_i - ib, block);
                A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);
                A2ASpatialSum<SpinColourVector_v> spatial_sum;
                spatial_sum.Allocate(Nii, Njj, grid);
                spatial_sum.PackLeftConj(left, ib, Nii);
                spatial_sum.PackRight(phaseGammaRight, 0, Njj);
                Eigen::Tensor<ComplexD, 3> block_result(nt, Nii, Njj);
                spatial_sum.Sum(block_result);
                for (int t=0;t<nt;t++) for (int ii=0;ii<Nii;ii++) for (int jj=0;jj<Njj;jj++)
                    mf(0, 0, t, ii, jj) = block_result(t, ii, jj);
#ifdef HADRONS_A2AM_PARALLEL_IO
                grid->Barrier();
                if (grid->ThisRank() == 0) {
#endif
                { A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j); io.saveBlock(mf, 0, 0, ib, jb); }
#ifdef HADRONS_A2AM_PARALLEL_IO
                } grid->Barrier();
#endif
            } // ib
        } // g, m
    } // jb
}
#endif // 0

#endif // Hadrons_MContraction_A2ANewMesonField_hpp_
