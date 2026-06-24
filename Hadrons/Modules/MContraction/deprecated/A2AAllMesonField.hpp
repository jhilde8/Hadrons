/*
 * A2AAllMesonField.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2024
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
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
#ifndef Hadrons_MContraction_A2AAllMesonField_hpp_
#define Hadrons_MContraction_A2AAllMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *  All-to-all meson field - spin-split BLAS, single GEMM per momentum block.
 *
 *  Pack left/right as [N*Ns, nxyz*Nc] instead of [N, nxyz*Nsc].  The GEMM
 *  output is a SpinMatrix [nt, Nii*Ns, Njj*Ns]; gamma is applied afterward as
 *  a cheap host-side trace, so one GEMM serves all gamma matrices at once.
 *  Phase trick: pack right at zero momentum, apply phase in-place before GEMM.
 *
 *  Loop order (jb, ib, m):
 *    PackLeftConjSpin, PackRightSpin  - once per (jb, ib)
 *    ApplyPhaseRight + Sum + Restore  - once per (jb, ib, m)
 *    gamma trace (host, cheap)        - all g per (jb, ib, m)
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

// Extends A2ASpatialSum with spin-split packing (N*Ns rows instead of N rows),
// enabling a single GEMM per momentum whose output carries all spin indices.
// This code lives here rather than in Grid because it is only used by this
// deprecated module; A2ASpatialSum itself remains general-purpose.
template<class vobj>
class A2ASpatialSumSpin : public A2ASpatialSum<vobj>
{
public:
    using typename A2ASpatialSum<vobj>::scalar;
    using typename A2ASpatialSum<vobj>::sobj;

    // Like Allocate, but splits each SpinColourVector into Ns spin rows.
    // N_i = _Nii * Ns, N_j = _Njj * Ns, Nsc = Nc = 3.
    // GEMM output [nt, Nii*Ns, Njj*Ns] enables post-GEMM gamma trace over all gammas.
    void AllocateSpin(int _Nii, int _Njj, GridBase *_grid)
    {
        const int Ns_qcd = 4;
        this->grid  = _grid;
        this->N_i   = _Nii * Ns_qcd;
        this->N_j   = _Njj * Ns_qcd;
        Coordinate ldims = this->grid->LocalDimensions();
        this->nt    = ldims[this->grid->Nd() - 1];
        this->nxyz  = this->grid->lSites() / this->nt;
        int Nsc_full = sizeof(sobj) / sizeof(scalar); // = 12 for SpinColourVector
        this->Nsc   = Nsc_full / Ns_qcd;              // = 3 = Nc

        this->W_buf.resize(this->nt * this->N_i * this->nxyz * this->Nsc);
        this->LR_buf.resize(this->nt * this->N_j * this->nxyz * this->Nsc);
        this->EMF_buf.resize(this->nt * this->N_j * this->N_i);

        this->W_ptrs.resize(this->nt);
        this->LR_ptrs.resize(this->nt);
        this->EMF_ptrs.resize(this->nt);
        scalar *Wh   = &this->W_buf[0];
        scalar *LRh  = &this->LR_buf[0];
        scalar *EMFh = &this->EMF_buf[0];
        int lN_i = this->N_i, lN_j = this->N_j, lnxyz = this->nxyz, lNsc = this->Nsc;
        for (int t = 0; t < this->nt; t++) {
            acceleratorPut(this->W_ptrs[t],   Wh   + t * lN_i * lnxyz * lNsc);
            acceleratorPut(this->LR_ptrs[t],  LRh  + t * lN_j * lnxyz * lNsc);
            acceleratorPut(this->EMF_ptrs[t], EMFh + t * lN_j * lN_i);
        }
    }

    void PackLeftConjSpin(const std::vector<Lattice<vobj>> &left, int start = 0, int count = -1)
    {
        const int Ns_qcd = 4;
        if (count < 0) count = (int)left.size();
        GRID_ASSERT(start + count <= (int)left.size());
        GRID_ASSERT(count * Ns_qcd == this->N_i);
        PackVectorsSpin<true>(left, &this->W_buf[0], count, start);
    }

    void PackRightSpin(const std::vector<Lattice<vobj>> &right, int start = 0, int count = -1)
    {
        const int Ns_qcd = 4;
        if (count < 0) count = (int)right.size();
        GRID_ASSERT(start + count <= (int)right.size());
        GRID_ASSERT(count * Ns_qcd == this->N_j);
        PackVectorsSpin<false>(right, &this->LR_buf[0], count, start);
    }

    // Pack N spin-colour vectors from vecs[start..start+N-1] into buf[nt][N*Ns][nxyz*Nc].
    // Mode n, spin s1 -> row (n*Ns + s1); color elements indexed by l_xyz*Nc + c.
    // DoConj=true conjugates each element during extraction (for PackLeftConjSpin).
    template<bool DoConj = false>
    void PackVectorsSpin(const std::vector<Lattice<vobj>> &vecs, scalar *buf, int N, int start = 0)
    {
        const int Ns_qcd = 4;
        int nd     = this->grid->_ndimension;
        int osites = this->grid->oSites();
        int Nsimd  = vobj::Nsimd();
        int lN_tot = N * Ns_qcd;
        int lNc    = this->Nsc;
        int lnxyz  = this->nxyz;
        Coordinate rdimensions = this->grid->_rdimensions;
        Coordinate ldims       = this->grid->LocalDimensions();
        Coordinate simd        = this->grid->_simd_layout;

        for (int n = 0; n < N; n++) {
            autoView(src_v, vecs[start + n], AcceleratorRead);
            accelerator_for(sf, osites, Nsimd, {
#ifdef GRID_SIMT
            {
                int lane = acceleratorSIMTlane(Nsimd);
#else
                for (int lane = 0; lane < Nsimd; lane++) {
#endif
                Coordinate icoor(nd), ocoor(nd), lcoor(nd);
                Lexicographic::CoorFromIndex(icoor, lane, simd);
                Lexicographic::CoorFromIndex(ocoor, sf, rdimensions);
                for (int d = 0; d < nd; d++)
                    lcoor[d] = rdimensions[d] * icoor[d] + ocoor[d];

                int     l_t = lcoor[nd - 1];
                Coordinate xyz_coor = lcoor;
                xyz_coor[nd - 1] = 0;
                int64_t l_xyz;
                Lexicographic::IndexFromCoor(xyz_coor, l_xyz, ldims);

                sobj    data   = extractLane(lane, src_v[sf]);
                if constexpr (DoConj) data = conjugate(data);
                scalar *data_s = (scalar *)&data;

                for (int s1 = 0; s1 < Ns_qcd; s1++) {
                    int64_t row  = (int64_t)n * Ns_qcd + s1;
                    int64_t base = (int64_t)l_t * lN_tot * lnxyz * lNc
                                 + row           * lnxyz  * lNc
                                 + l_xyz         * lNc;
                    for (int c = 0; c < lNc; c++)
                        buf[base + c] = data_s[s1 * lNc + c];
                }
            }
            });
        }
    }
};

class A2AAllMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AAllMesonFieldPar,
                                    int,                      block,
                                    std::string,              left,
                                    std::string,              right,
                                    std::string,              output,
                                    std::string,              gammas,
                                    std::vector<std::string>, mom);
};

class A2AAllMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AAllMesonFieldMetadata,
                                    std::vector<RealF>, momentum,
                                    Gamma::Algebra,     gamma);
};

template <typename FImpl>
class TA2AAllMesonField : public Module<A2AAllMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TA2AAllMesonField(const std::string name);
    virtual ~TA2AAllMesonField(void){};
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

MODULE_REGISTER(A2AAllMesonField, ARG(TA2AAllMesonField<FIMPL>), MContraction);

/******************************************************************************
 *                  TA2AAllMesonField implementation                          *
 ******************************************************************************/
template <typename FImpl>
TA2AAllMesonField<FImpl>::TA2AAllMesonField(const std::string name)
: Module<A2AAllMesonFieldPar>(name)
, momphName_(name + "_momph")
{}

template <typename FImpl>
std::vector<std::string> TA2AAllMesonField<FImpl>::getInput(void)
{
    return {par().left, par().right};
}

template <typename FImpl>
std::vector<std::string> TA2AAllMesonField<FImpl>::getOutput(void)
{
    return {};
}

template <typename FImpl>
void TA2AAllMesonField<FImpl>::setup(void)
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
void TA2AAllMesonField<FImpl>::execute(void)
{
    typedef typename FImpl::SiteSpinor      vobj;
    typedef typename vobj::vector_type      vector_type;
    typedef iSpinColourVector<vector_type>  SpinColourVector_v;
    typedef typename SpinColourVector_v::scalar_type scalar_t;

    const int Ns = 4; // QCD spinor dimension

    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = envGetGrid(FermionField);

    int nt     = env().getDim().back();
    int N_i    = left.size();
    int N_j    = right.size();
    int ngamma = gamma_.size();
    int nmom   = mom_.size();
    int block  = par().block;

    LOG(Message) << "Computing all-to-all meson fields (spin-split BLAS)" << std::endl;
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
        A2AAllMesonFieldMetadata md;
        for (auto pmu: mom_[m])
            md.momentum.push_back(pmu);
        md.gamma = gamma_[g];
        return md;
    };

    // Output buffer: one (nt, Nii, Njj) block at a time (reused per gamma).
    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

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
        A2AAllMesonFieldMetadata md = metadataFn(m, g);
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

    // Pre-pack flat phase arrays for in-place LR_buf multiplication.
    std::vector<deviceVector<scalar_t>> ph_flat(nmom);
    for (int m = 0; m < nmom; m++)
        A2ASpatialSum<SpinColourVector_v>::PackPhase(grid, ph[m], ph_flat[m]);

    // Loop order (jb, ib, m):
    //   PackLeftConjSpin, PackRightSpin - once per (jb, ib)
    //   ApplyPhaseRight + Sum + Restore - once per (jb, ib, m)
    //   gamma trace (host, all g)       - once per (jb, ib, m)

    for (int jb = 0; jb < N_j; jb += block)
    {
        int Njj = std::min(N_j - jb, block);

        for (int ib = 0; ib < N_i; ib += block)
        {
            int Nii = std::min(N_i - ib, block);

            A2ASpatialSumSpin<SpinColourVector_v> spatial_sum;
            spatial_sum.AllocateSpin(Nii, Njj, grid);
            spatial_sum.PackLeftConjSpin(left,  ib, Nii);
            spatial_sum.PackRightSpin(right, jb, Njj);

            // Reused across all momenta; size is fixed for this (ib, jb) block.
            Eigen::Tensor<ComplexD, 3> spin_result(nt, Nii * Ns, Njj * Ns);

            for (int m = 0; m < nmom; m++)
            {
                spatial_sum.ApplyPhaseRight(ph_flat[m]);
                spatial_sum.Sum(spin_result);
                spatial_sum.RestorePhaseRight(ph_flat[m]);

                for (int g = 0; g < ngamma; g++)
                {
                    A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);
                    for (int t  = 0; t  < nt;  t++)
                    for (int ii = 0; ii < Nii; ii++)
                    for (int jj = 0; jj < Njj; jj++)
                    {
                        // MF(ii,jj,t) = sum_{sl,sr} G[sl,sr] * spin_result[ii*Ns+sl, jj*Ns+sr]
                        //             = trace(G * M) where M[sr,sl] = spin_result[ii*Ns+sl, jj*Ns+sr].
                        // Transposed fill: trace(G*M) = sum_{a,b} G[a,b]*M[b,a],
                        // so M[b,a] = spin_result[..+a, ..+b] requires spinMat(sr,sl) = spin_result(sl,sr).
                        // All Ns*Ns entries are written below, so no Zero() initialisation needed.
                        SpinMatrixD spinMat;
                        for (int sl = 0; sl < Ns; sl++)
                        for (int sr = 0; sr < Ns; sr++)
                            spinMat()(sr, sl)() = spin_result(t, ii*Ns + sl, jj*Ns + sr);
                        mf(0, 0, t, ii, jj) = TensorRemove(trace(Gamma(gamma_[g]) * spinMat));
                    }

                    std::string ioname   = ionameFn(m, g);
                    std::string filename = filenameFn(m, g);

                    LOG(Message) << "AMF block i=" << ib << " j=" << jb
                                 << " m=" << m << " g=" << gamma_[g] << std::endl;

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
                } // g
            } // m
        } // ib
    } // jb
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AAllMesonField_hpp_
