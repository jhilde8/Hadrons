/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2AChromoMagneticOperatorField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_
#define Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *          All-to-all chromo-magnetic operator field creation (GPU)          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AChromoMagneticOperatorFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldPar,
                                    int,         cacheBlock,
                                    std::string, parities,
                                    std::string, left,
                                    std::string, right,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, ifOrthogs);
};

class A2AChromoMagneticOperatorFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldMetadata,
                                    std::string, meta);
};

template <typename GImpl, typename FImpl>
class TA2AChromoMagneticOperatorField : public Module<A2AChromoMagneticOperatorFieldPar>
{
public:
  typedef typename GImpl::GaugeLinkField GaugeMat;
  typedef typename GImpl::SiteGaugeLink  SiteGaugeLink;
  FERM_TYPE_ALIASES(FImpl,);
  typedef typename FImpl::SiteSpinor     vobj;
  typedef typename vobj::scalar_object   sobj;
  typedef typename vobj::scalar_type     scalar_type;
  typedef typename vobj::vector_type     vector_type;
public:
    // constructor
    TA2AChromoMagneticOperatorField(const std::string name);
    // destructor
    virtual ~TA2AChromoMagneticOperatorField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::vector<int> parities_;
  std::vector<int> ifOrthogs_;
};

MODULE_REGISTER(A2AChromoMagneticOperatorField, ARG(TA2AChromoMagneticOperatorField<GIMPL,FIMPL>), MContraction);

/******************************************************************************
*               TA2AChromoMagneticOperatorField implementation                *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
TA2AChromoMagneticOperatorField<GImpl,FImpl>::TA2AChromoMagneticOperatorField(const std::string name)
: Module<A2AChromoMagneticOperatorFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
std::vector<std::string> TA2AChromoMagneticOperatorField<GImpl,FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().gauge};
  return in;
}

template <typename GImpl, typename FImpl>
std::vector<std::string> TA2AChromoMagneticOperatorField<GImpl,FImpl>::getOutput(void)
{
  std::vector<std::string> out = {};
  return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2AChromoMagneticOperatorField<GImpl,FImpl>::setup(void)
{
  parities_  = strToVec<int>(par().parities);
  ifOrthogs_ = strToVec<int>(par().ifOrthogs);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename GImpl, typename FImpl>
void TA2AChromoMagneticOperatorField<GImpl,FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iColourMatrix<vector_type>     ColourMatrix_v;
  typedef iSpinMatrix<vector_type>       SpinMatrix_v;
  typedef iSinglet<vector_type>          Scalar_v;
  typedef iSinglet<scalar_type>          Scalar_s;

  auto &left    = envGet(std::vector<FermionField>, par().left);
  auto &right   = envGet(std::vector<FermionField>, par().right);
  const auto &U = envGet(GaugeField, par().gauge);

  GridBase *grid = left[0].Grid();

  int orthogdim = 3;
  int rd     = grid->_rdimensions[orthogdim];
  int ld     = grid->_ldimensions[orthogdim];
  int Nd     = grid->_ndimension;
  int Nsimd  = grid->Nsimd();
  int e1     = grid->_slice_nblock[orthogdim];
  int e2     = grid->_slice_block [orthogdim];
  int stride = grid->_slice_stride[orthogdim];

  LOG(Message) << "Computing all-to-all ChromoMagnetic operator fields (GPU)" << std::endl;
  LOG(Message) << "R dimension: "  << rd     << std::endl;
  LOG(Message) << "Slice nblock: " << e1     << std::endl;
  LOG(Message) << "Slice block: "  << e2     << std::endl;
  LOG(Message) << "Slice stride: " << stride << std::endl;

  int nt         = env().getDim().back();
  int N_i        = left.size();
  int N_j        = right.size();
  int cacheBlock = par().cacheBlock;
  Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(nt*N_i*N_j);

  LOG(Message) << "Left: '"        << par().left  << "' Right: '"
               << par().right      << "'"          << std::endl;
  LOG(Message) << "Gauge field: '" << par().gauge  << "'"          << std::endl;
  LOG(Message) << "Parities:"      << std::endl;
  for (auto &p: parities_)
    LOG(Message) << "  " << p << std::endl;
  LOG(Message) << "CMO field size: " << nt << "*" << N_i << "*" << N_j
               << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) << std::endl;

  double t0, t1, t2, t3;

  std::vector<FermionField> loopRight(cacheBlock, grid);

  for (auto &ifOrthog: ifOrthogs_) {
    std::vector<GaugeMat>  G;
    Vector<Gamma::Algebra> Sigma;

    t0 = -usecond();
    if (ifOrthog == 0)
      Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContraction0(G, Sigma, U);
    else
      Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContraction1(G, Sigma, U);
    t0 += usecond();
    LOG(Message) << "Field strength constructed for ifOrthog=" << ifOrthog
                 << "; t_fs=" << t0 << " us" << std::endl;

    for (auto &parity: parities_) {
      LOG(Message) << "Starting calculation with ifOrthog=" << ifOrthog
                   << " parity=" << parity << std::endl;
      A2AMatrixSet<HADRONS_A2AM_IO_TYPE> cmf(mBuf.data(), 1, 1, nt, N_i, N_j);
      A2ASpatialSum<SpinColourVector_v> spatial_sum;

      t1 = -usecond();
      LOG(Message) << "Making CMF" << std::endl;

      for (unsigned int j = 0; j < N_j; j += cacheBlock) {
        int Njj = MIN(N_j-j, cacheBlock);

        t2 = -usecond();
        for (int jj = 0; jj < Njj; jj++)
          Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContractRight(
              loopRight[jj], G, Sigma, right[j+jj], parity);
        t2 += usecond();
        LOG(Message) << "loopRight packed for j-block " << j/cacheBlock
                     << "; ifOrthog=" << ifOrthog << " parity=" << parity
                     << "; t_LR=" << t2 << " us" << std::endl;

        t3 = -usecond();
        for (unsigned int i = 0; i < N_i; i += cacheBlock) {
          int Nii = MIN(N_i-i, cacheBlock);

          spatial_sum.Allocate(Nii, Njj, grid);
          spatial_sum.PackLeftConj(left, i, Nii);
          spatial_sum.PackRight(loopRight, 0, Njj);

          Eigen::Tensor<ComplexD,3> cmfBlock(nt, Nii, Njj);
          cmfBlock.setZero();
          spatial_sum.Sum(cmfBlock);

          for (int t  = 0; t  < nt;  t++)
          for (int ii = 0; ii < Nii; ii++)
          for (int jj = 0; jj < Njj; jj++)
            cmf(0,0,t,i+ii,j+jj) = cmfBlock(t,ii,jj);

          LOG(Message) << "CMF made for i-block " << i/cacheBlock
                       << "; j-block "            << j/cacheBlock
                       << "; ifOrthog="           << ifOrthog
                       << " parity="              << parity << std::endl;
        }// i
        t3 += usecond();
        LOG(Message) << "CMF made for all i-blocks in j-block " << j/cacheBlock
                     << "; t_block=" << t3 << " us" << std::endl;
      }// j
      t1 += usecond();
      LOG(Message) << "CMF made; t_cmf=" << t1 << " us" << std::endl;

#ifdef HADRONS_A2AM_PARALLEL_IO
      grid->Barrier();
      if (grid->ThisRank() == 0) {
#endif
      std::string ioname = "parity" + std::to_string(parity);
      if (ifOrthog == 1)
        ioname = ioname + "_GijSij";
      else
        ioname = ioname + "_GitSit";
      std::string filename = par().output + "." + std::to_string(vm().getTrajectory())
                           + "/" + ioname + ".h5";
      LOG(Message) << "Writing block to " << filename << std::endl;
      makeFileDir(filename, grid);
      A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename, ioname, nt, N_i, N_j);
      A2AChromoMagneticOperatorFieldMetadata md;
      md.meta = ioname;
      io.initFile(md, MAX(N_i,N_j));
      io.saveBlock(cmf, 0, 0, 0, 0);
#ifdef HADRONS_A2AM_PARALLEL_IO
      }
      grid->Barrier();
#endif
    }// parity
  }// ifOrthog
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_
