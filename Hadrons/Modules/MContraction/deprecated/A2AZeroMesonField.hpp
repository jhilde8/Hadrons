/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: Hadrons/Modules/MContraction/A2AZeroMesonField.hpp

    Copyright (C) 2015-2024

*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AZeroMesonField_hpp_
#define Hadrons_MContraction_A2AZeroMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

class A2AZeroMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AZeroMesonFieldPar,
                                    int,         block,
                                    int,         a2aBlock,
                                    std::string, left,
                                    std::string, right,
                                    std::string, output,
                                    std::string, gammas);
};

class A2AZeroMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AZeroMesonFieldMetadata,
                                    Gamma::Algebra, gamma);
};

template <typename FImpl>
class TA2AZeroMesonField : public Module<A2AZeroMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    TA2AZeroMesonField(const std::string name);
    virtual ~TA2AZeroMesonField(void){};
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual void setup(void);
    virtual void execute(void);
private:
    std::vector<Gamma::Algebra> gamma_;

    // Zero-momentum meson field via batched GEMM (no MomentumProject).
    // Moved here from A2Autils so A2ASpatialSum.h stays free of ZMF-specific code.
    static void ZeroMesonField(Eigen::Tensor<ComplexD, 3> &mat,
                               const std::vector<FermionField> &lhs_wi,
                               int lhs_start, int lhs_count,
                               const std::vector<FermionField> &rhs_vj,
                               int rhs_start, int rhs_count,
                               Gamma::Algebra gamma,
                               int a2aBlock)
    {
        typedef iSpinColourVector<typename FImpl::SiteSpinor::vector_type> SpinColourVector_v;

        GridBase *grid      = lhs_wi[lhs_start].Grid();
        int       Ni        = lhs_count;
        int       Nj        = rhs_count;
        int       Nsimd     = grid->Nsimd();
        uint64_t  oSites    = grid->oSites();
        int       nt_global = (int)mat.dimension(0);

        double t_gamma = 0;
        double t_gemm  = 0;

        for (int jo = 0; jo < Nj; jo += a2aBlock) {
            int Njj = std::min(Nj - jo, a2aBlock);

            std::vector<FermionField> gammaRight(Njj, grid);
            {
                Gamma::Algebra ga = gamma;
                for (int jj = 0; jj < Njj; jj++) {
                    t_gamma -= usecond();
                    autoView(outv, gammaRight[jj],                AcceleratorWrite);
                    autoView(inv,  rhs_vj[rhs_start + jo + jj],  AcceleratorRead);
                    accelerator_for(ss, oSites, (size_t)Nsimd, {
                        coalescedWrite(outv[ss], Gamma(ga) * inv(ss));
                    });
                    t_gamma += usecond();
                }
            }

            for (int io = 0; io < Ni; io += a2aBlock) {
                int Nii = std::min(Ni - io, a2aBlock);

                t_gemm -= usecond();
                A2ASpatialSum<SpinColourVector_v> spatial_sum;
                spatial_sum.Allocate    (Nii, Njj, grid);
                spatial_sum.PackLeftConj(lhs_wi,     lhs_start + io, Nii);
                spatial_sum.PackRight   (gammaRight, 0,              Njj);

                Eigen::Tensor<ComplexD, 3> mfChunk(nt_global, Nii, Njj);
                spatial_sum.Sum(mfChunk);
                t_gemm += usecond();

                for (int t  = 0; t  < nt_global; t++)
                for (int ii = 0; ii < Nii;       ii++)
                for (int jj = 0; jj < Njj;       jj++)
                    mat(t, io + ii, jo + jj) = mfChunk(t, ii, jj);
            }
        }

        std::cout << GridLogMessage << " ZeroMesonField t_gamma " << t_gamma/1.e6 << "s" << std::endl;
        std::cout << GridLogMessage << " ZeroMesonField t_gemm  " << t_gemm/1.e6  << "s" << std::endl;
    }
};

MODULE_REGISTER(A2AZeroMesonField, ARG(TA2AZeroMesonField<FIMPL>), MContraction);

/******************************************************************************
 *                 TA2AZeroMesonField implementation                          *
 ******************************************************************************/
template <typename FImpl>
TA2AZeroMesonField<FImpl>::TA2AZeroMesonField(const std::string name)
: Module<A2AZeroMesonFieldPar>(name)
{}

template <typename FImpl>
std::vector<std::string> TA2AZeroMesonField<FImpl>::getInput(void)
{
    return {par().left, par().right};
}

template <typename FImpl>
std::vector<std::string> TA2AZeroMesonField<FImpl>::getOutput(void)
{
    return {};
}

template <typename FImpl>
void TA2AZeroMesonField<FImpl>::setup(void)
{
    gamma_.clear();
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
}

template <typename FImpl>
void TA2AZeroMesonField<FImpl>::execute(void)
{
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);

    GridBase *grid = left[0].Grid();

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int block    = par().block;
    int a2aBlock = par().a2aBlock;

    LOG(Message) << "Computing all-to-all zero-momentum meson fields" << std::endl;
    LOG(Message) << "Left: '" << par().left << "' Right: '" << par().right << "'" << std::endl;
    LOG(Message) << "Spin bilinears:" << std::endl;
    for (auto &g: gamma_)
        LOG(Message) << "  " << g << std::endl;
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE))
                 << "/bilinear)" << std::endl;

    // Output buffer for one (Nii x Njj x nt) block at a time.
    Vector<HADRONS_A2AM_IO_TYPE> mBuf;
    mBuf.resize(nt * block * block);

    // Create output directory once before the gamma loop.
    std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());
    {
        std::string dummy = dirBase + "/mkdir.h5";
        makeFileDir(dummy, grid);
    }

    for (auto &gamma: gamma_)
    {
        std::stringstream ionameSS;
        ionameSS << gamma << "_0_0_0";
        std::string ioname   = ionameSS.str();
        std::string filename = dirBase + "/" + ioname + ".h5";

        A2AZeroMesonFieldMetadata md;
        md.gamma = gamma;

        // Initialise HDF5 file for this gamma before the block loops.
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

        // Outer j-block loop.
        for (int jb = 0; jb < N_j; jb += block)
        {
            int Njj = std::min(N_j - jb, block);

            // Outer i-block loop.
            for (int ib = 0; ib < N_i; ib += block)
            {
                int Nii = std::min(N_i - ib, block);

                // mf is a non-owning view over the block-sized buffer.
                A2AMatrixSet<HADRONS_A2AM_IO_TYPE> mf(mBuf.data(), 1, 1, nt, Nii, Njj);

                double t_block = -usecond();

                Eigen::Tensor<ComplexD, 3> mfResult(nt, Nii, Njj);
                ZeroMesonField(
                    mfResult,
                    left,  ib, Nii,
                    right, jb, Njj,
                    gamma, a2aBlock);

                for (int t  = 0; t  < nt;  t++)
                for (int ii = 0; ii < Nii; ii++)
                for (int jj = 0; jj < Njj; jj++)
                    mf(0, 0, t, ii, jj) = mfResult(t, ii, jj);

                t_block += usecond();
                LOG(Message) << "MF block i=" << ib << " j=" << jb
                             << " gamma=" << gamma
                             << "; t=" << t_block/1.e6 << " s" << std::endl;

                // Write this (ib, jb) block to HDF5.
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
            } // ib
        } // jb
    } // gamma
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AZeroMesonField_hpp_
