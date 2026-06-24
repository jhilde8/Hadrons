/*
 * A2AMomMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
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
#ifndef Hadrons_MContraction_A2AMomMesonField_hpp_
#define Hadrons_MContraction_A2AMomMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *   All-to-all meson field: momentum-factored two-phase custom kernel path.
 *   Phase 1 (SpinColorTrace): SC inner product with on-the-fly gamma, no phase.
 *   Phase 2 (MomMesonField):  spatial sum against all momenta simultaneously.
 *   Eliminates the per-momentum SC contraction redundancy of A2AMesonField
 *   and the gammaRight temporary.
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AMomMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMomMesonFieldPar,
                                    int,                      block,
                                    std::string,              left,
                                    std::string,              right,
                                    std::string,              output,
                                    std::string,              gammas,
                                    std::vector<std::string>, mom);
};

class A2AMomMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AMomMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma);
};

template <typename FImpl>
class TA2AMomMesonField : public Module<A2AMomMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TA2AMomMesonField(const std::string name);
    virtual ~TA2AMomMesonField(void){};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
private:
    bool                           hasPhase_{false};
    std::string                    momphName_;
    std::vector<Gamma::Algebra>    gamma_;
    std::vector<std::vector<Real>> mom_;
};

MODULE_REGISTER(A2AMomMesonField, ARG(TA2AMomMesonField<FIMPL>), MContraction);

/******************************************************************************
 *                  TA2AMomMesonField implementation                          *
 ******************************************************************************/
template <typename FImpl>
TA2AMomMesonField<FImpl>::TA2AMomMesonField(const std::string name)
: Module<A2AMomMesonFieldPar>(name)
, momphName_(name + "_momph")
{}

template <typename FImpl>
std::vector<std::string> TA2AMomMesonField<FImpl>::getInput(void)
{
    return {par().left, par().right};
}

template <typename FImpl>
std::vector<std::string> TA2AMomMesonField<FImpl>::getOutput(void)
{
    return {};
}

template <typename FImpl>
void TA2AMomMesonField<FImpl>::setup(void)
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
void TA2AMomMesonField<FImpl>::execute(void)
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

    int nd       = grid->Nd();
    int nt_local = grid->LocalDimensions()[nd - 1];
    int nxyz     = grid->lSites() / nt_local;
    int Nsc      = sizeof(typename SpinColourVector_v::scalar_object) / sizeof(scalar_t);

    {
        int64_t q_bytes = (int64_t)block * block * nt_local * nxyz * sizeof(scalar_t);
        if (q_bytes >= (int64_t)60 * 1024 * 1024 * 1024)
            HADRONS_ERROR(Size, "q_buf would exceed 60 GB; reduce par().block");
    }

    LOG(Message) << "Computing all-to-all meson fields (momentum-factored kernel)"
                 << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'"
                 << std::endl;
    LOG(Message) << "Momenta:" << std::endl;
    for (auto &p: mom_)
        LOG(Message) << "  " << p << std::endl;
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
        LOG(Message) << "  " << g << std::endl;
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize "
                 << sizeString(nt * N_i * N_j * sizeof(HADRONS_A2AM_IO_TYPE))
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
                ph[j] = ph[j] + (mom_[j][mu] / env().getDim(mu)) * coor;
            }
            ph[j] = exp((Real)(2 * M_PI) * i * ph[j]);
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
        A2AMomMesonFieldMetadata md;
        for (auto pmu: mom_[m])
            md.momentum.push_back(pmu);
        md.gamma = gamma_[g];
        return md;
    };

    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

    std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());
    {
        std::string dummy = dirBase + "/mkdir.h5";
        makeFileDir(dummy, grid);
    }

    for (int m = 0; m < nmom; m++)
    for (int g = 0; g < ngamma; g++)
    {
        std::string filename = filenameFn(m, g);
        std::string ioname   = ionameFn(m, g);
        A2AMomMesonFieldMetadata md = metadataFn(m, g);
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

    // Pack all nmom phases into one flat buffer [m*nxyz + x], done once.
    startTimer("Pack phases");
    deviceVector<scalar_t> phase_all;
    A2Autils<FImpl>::PackAllPhases(grid, ph, phase_all);
    stopTimer("Pack phases");

    // Loop order (jb, g, ib):
    //   PackLeftConj   - once per (jb, g, ib)
    //   PackRight      - once per (jb, g, ib); packs right directly, no gamma applied
    //   SpinColorTrace - once per (jb, g, ib); gamma applied on-the-fly, all sites
    //   MomMesonField  - once per (jb, g, ib); all nmom momenta simultaneously
    //   IO             - nmom writes per (jb, g, ib)

    for (int jb = 0; jb < N_j; jb += block)
    {
        int Njj = std::min(N_j - jb, block);

        for (int g = 0; g < ngamma; g++)
        {
            for (int ib = 0; ib < N_i; ib += block)
            {
                int Nii = std::min(N_i - ib, block);

                A2ASpatialSum<SpinColourVector_v> spatial_sum;
                spatial_sum.Allocate(Nii, Njj, grid);

                startTimer("Pack vectors");
                spatial_sum.PackLeftConj(left, ib, Nii);
                spatial_sum.PackRight(right, jb, Njj);
                stopTimer("Pack vectors");

                deviceVector<scalar_t> q_buf;
                startTimer("SpinColorTrace");
                A2Autils<FImpl>::SpinColorTrace(
                    q_buf,
                    spatial_sum.W_buf,
                    spatial_sum.LR_buf,
                    gamma_[g],
                    Nii, Njj, nt_local, nxyz, Nsc);
                stopTimer("SpinColorTrace");

                Eigen::Tensor<ComplexD, 4> block_result(nt, Nii, Njj, nmom);
                startTimer("MomMesonField");
                A2Autils<FImpl>::MomMesonField(
                    block_result,
                    q_buf,
                    phase_all,
                    Nii, Njj, nt_local, nxyz, nmom,
                    grid);
                stopTimer("MomMesonField");

                startTimer("IO");
                for (int m = 0; m < nmom; m++)
                {
                    A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);
                    for (int t  = 0; t  < nt;  t++)
                    for (int ii = 0; ii < Nii; ii++)
                    for (int jj = 0; jj < Njj; jj++)
                        mf(0, 0, t, ii, jj) = block_result(t, ii, jj, m);

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
                } // m
                stopTimer("IO");
            } // ib
        } // g
    } // jb
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AMomMesonField_hpp_
