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
                                    int,         cacheBlock,
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
    int block      = par().block;
    int cacheBlock = par().cacheBlock;

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

                // cacheBlock sub-loops: one ZeroMesonField call per sub-block.
                for (int sub_jj = 0; sub_jj < Njj; sub_jj += cacheBlock)
                {
                    int Njjj = std::min(Njj - sub_jj, cacheBlock);
                    for (int sub_ii = 0; sub_ii < Nii; sub_ii += cacheBlock)
                    {
                        int Niii = std::min(Nii - sub_ii, cacheBlock);

                        Eigen::Tensor<ComplexD, 3> mfChunk(nt, Niii, Njjj);
                        A2Autils<FImpl>::ZeroMesonField(
                            mfChunk,
                            left,  ib + sub_ii, Niii,
                            right, jb + sub_jj, Njjj,
                            gamma);

                        for (int t   = 0; t   < nt;   t++)
                        for (int ii  = 0; ii  < Niii; ii++)
                        for (int jjj = 0; jjj < Njjj; jjj++)
                            mf(0, 0, t, sub_ii + ii, sub_jj + jjj) = mfChunk(t, ii, jjj);
                    }
                }

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
