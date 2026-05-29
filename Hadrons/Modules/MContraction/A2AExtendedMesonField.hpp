/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MContraction/A2AExtendedMesonField.hpp

Copyright (C) 2015-2019

Author: Masaaki Tomii <masaaki.tomii@uconn.edu>
*************************************************************************************/
/*  END LEGAL */
#ifndef Hadrons_MContraction_A2AExtendedMesonField_hpp_
#define Hadrons_MContraction_A2AExtendedMesonField_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/A2AMatrix.hpp>
#include <Grid/qcd/utils/A2Autils.h>

BEGIN_HADRONS_NAMESPACE
//   _
//  / \
//  \ /
// --x--
//

/******************************************************************************
 *                All-to-all extended meson field creation                    *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AExtendedMesonFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldPar,
                                    int, cacheBlock,
				    std::string, types,
                                    std::string, left,
                                    std::string, right,
				    std::string, loop_vw1,
				    std::string, loop_vw2,
                                    std::string, output,
                                    std::string, gammas1,
				    std::string, gammas2);
};

class A2AExtendedMesonFieldMetadata: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AExtendedMesonFieldMetadata,
                                    std::string, gamma1,
				    std::string, gamma2);
};

template <typename FImpl>
class TA2AExtendedMesonField : public Module<A2AExtendedMesonFieldPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
  //typedef typename FImpl::FermionField FermionField;
  typedef typename FImpl::SiteSpinor vobj;
  typedef typename vobj::scalar_object sobj;
  typedef typename vobj::scalar_type scalar_type;
  typedef typename vobj::vector_type vector_type;
public:
    // constructor
    TA2AExtendedMesonField(const std::string name);
    // destructor
    virtual ~TA2AExtendedMesonField(void){};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
  std::vector<std::vector<Gamma::Algebra> >       gamma1_;
  std::vector<std::vector<Gamma::Algebra> >       gamma2_;
  std::vector<std::string> nameg1_;
  std::vector<std::string> nameg2_;
  std::vector<int> types_;
};

MODULE_REGISTER(A2AExtendedMesonField, TA2AExtendedMesonField<FIMPL>, MContraction);

/******************************************************************************
*               TA2AExtendedMesonField implementation                         *
******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TA2AExtendedMesonField<FImpl>::TA2AExtendedMesonField(const std::string name)
: Module<A2AExtendedMesonFieldPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
// Modification needed??
template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonField<FImpl>::getInput(void)
{
  std::vector<std::string> in = {par().left, par().right, par().loop_vw1, par().loop_vw2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TA2AExtendedMesonField<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};

    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::setup(void)
{
    gamma1_.clear();
    gamma2_.clear();
    std::vector<std::string> tmp1 = strToVec<std::string>(par().gammas1);
    std::vector<std::string> tmp2 = strToVec<std::string>(par().gammas2);
    nameg1_ = tmp1;
    nameg2_ = tmp2;
    types_ = strToVec<int>(par().types);
    assert( tmp1.size() == tmp2.size() );
    for ( int ig = 0; ig < tmp1.size(); ++ig ) {
      std::vector<Gamma::Algebra> vec;
      vec.clear();
      if ( tmp1[ig] == "GammaMU" ) {
	vec = {
	       Gamma::Algebra::GammaX,
	       Gamma::Algebra::GammaY,
	       Gamma::Algebra::GammaZ,
	       Gamma::Algebra::GammaT
	};
      } else if ( tmp1[ig] == "GammaMUGamma5" ) {
	vec = {
	       Gamma::Algebra::GammaXGamma5,
	       Gamma::Algebra::GammaYGamma5,
	       Gamma::Algebra::GammaZGamma5,
	       Gamma::Algebra::GammaTGamma5
	};
      } else if ( tmp1[ig] == "SigmaMUNU" ) {
	vec = {
	       Gamma::Algebra::SigmaXY,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::SigmaXT,
	       Gamma::Algebra::SigmaYZ,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::SigmaZT
	};
      } else if ( tmp1[ig] == "SigmaMUNUGamma5" ) {
	vec = {
	       Gamma::Algebra::MinusSigmaZT,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::MinusSigmaYZ,
	       Gamma::Algebra::MinusSigmaXT,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::MinusSigmaXY
	};
      } else {
	vec = strToVec<Gamma::Algebra>(tmp1[ig]);
      }
      gamma1_.push_back(vec);

      vec.clear();
      if ( tmp2[ig] == "GammaMU" ) {
	vec = {
	       Gamma::Algebra::GammaX,
	       Gamma::Algebra::GammaY,
	       Gamma::Algebra::GammaZ,
	       Gamma::Algebra::GammaT
	};
      } else if ( tmp2[ig] == "GammaMUGamma5" ) {
	vec = {
	       Gamma::Algebra::GammaXGamma5,
	       Gamma::Algebra::GammaYGamma5,
	       Gamma::Algebra::GammaZGamma5,
	       Gamma::Algebra::GammaTGamma5
	};
      } else if ( tmp2[ig] == "SigmaMUNU" ) {
	vec = {
	       Gamma::Algebra::SigmaXY,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::SigmaXT,
	       Gamma::Algebra::SigmaYZ,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::SigmaZT
	};
      } else if ( tmp2[ig] == "SigmaMUNUGamma5" ) {
	vec = {
	       Gamma::Algebra::MinusSigmaZT,
	       Gamma::Algebra::SigmaYT,
	       Gamma::Algebra::MinusSigmaYZ,
	       Gamma::Algebra::MinusSigmaXT,
	       Gamma::Algebra::SigmaXZ,
	       Gamma::Algebra::MinusSigmaXY
	};
      } else {
	vec = strToVec<Gamma::Algebra>(tmp2[ig]);
      }
      gamma2_.push_back(vec);
    }
    auto &left  = envGet(std::vector<FermionField>, par().left);
    GridBase *grid = left[0].Grid();
    envCreateLat(PropagatorField, "propagatorLoop");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
  typedef iSpinColourMatrix<vector_type> SpinColourMatrix_v;
  typedef iSpinMatrix<vector_type> SpinMatrix_v;
  typedef iSinglet<vector_type> Scalar_v;
  typedef iSinglet<scalar_type> Scalar_s;
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &loop1 = envGet(std::vector<FermionField>, par().loop_vw1);
    auto &loop2 = envGet(std::vector<FermionField>, par().loop_vw2);

    GridBase *grid = left[0].Grid();

    int orthogdim = 3;
    int rd=grid->_rdimensions[orthogdim];//2
    int ld=grid->_ldimensions[orthogdim];
    int Nd=grid->_ndimension;
    int Nsimd=grid->Nsimd();
    int e1=    grid->_slice_nblock[orthogdim];//1
    int e2=    grid->_slice_block [orthogdim];//64 must be 4^3
    int stride=grid->_slice_stride[orthogdim];//128
    LOG(Message) << "Computing all-to-all EXTENDED meson fields" << std::endl;
    LOG(Message) << "R dimension: " << rd << std::endl;
    LOG(Message) << "Slice nblock: " << e1 << std::endl;
    LOG(Message) << "Slice block: " << e2 << std::endl;
    LOG(Message) << "Slice stride: " << stride << std::endl;

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    //int block      = par().block;
    int cacheBlock = par().cacheBlock;
    Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(nt*N_i*N_j);

  int ngamma = 0;
    for ( int ig = 0 ; ig < gamma1_.size() ; ++ig ) {
      ngamma += gamma1_[ig].size();
    }

    LOG(Message) << "Left: '" << par().left << "' Right: '"
		 << par().right << "'" << std::endl;
    LOG(Message) << "A2AVectors for loop: '" << par().loop_vw1
		 << "' and '" << par().loop_vw2 << "'" << std::endl;
    LOG(Message) << "Spin bilinears1:" << std::endl;
    for (auto &g: gamma1_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Spin bilinears2:" << std::endl;
    for (auto &g: gamma2_)
    {
        LOG(Message) << "  " << g << std::endl;
    }
    LOG(Message) << "Meson field size: " << nt << "*" << N_i << "*" << N_j 
                 << " (filesize " << sizeString(nt*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE)) 
                 << "/momentum/bilinear)" << std::endl;

    auto &loop = envGet(PropagatorField, "propagatorLoop");
    Grid::A2AExtendedMesonField<FImpl>::LoopPropagator(loop, loop1, loop2);
    LOG(Message) << "Quark loop calculated" << std::endl;

    std::vector<FermionField> leftv(N_i, grid);
    for (int i = 0; i < N_i; i++)
      Grid::A2AExtendedMesonField<FImpl>::PackLeftConjugated(leftv[i], left[i]);
    LOG(Message) << "Memory for left vectors allocated" << std::endl;

    for (int &type: types_){

      for (int ig = 0 ; ig < gamma1_.size() ; ++ig ){
		
	A2AMatrixSet<HADRONS_A2AM_IO_TYPE> emf(mBuf.data(),1,1,nt,N_i,N_j);
	
	Vector<Gamma::Algebra> gamma1(gamma1_[ig].begin(), gamma1_[ig].end());
	Vector<Gamma::Algebra> gamma2(gamma2_[ig].begin(), gamma2_[ig].end());

	PropagatorField tloop(grid);
	tloop = Zero();
	switch (type) {
	case 0: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType0(tloop, loop);                 break;
	case 1: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType1(tloop, loop, gamma1, gamma2); break;
	case 2: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType2(tloop, loop, gamma2);         break;
	case 3: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType3(tloop, loop, gamma1, gamma2); break;
	}
	LOG(Message) << "tloop contraction done for type " << type << std::endl;

	std::vector<FermionField> loopRight(N_j, grid);
	for (int j = 0; j < N_j; j++) {
	  switch (type) {
	  case 0: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType0(loopRight[j], tloop, right[j], gamma1, gamma2); break;
	  case 1: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType1(loopRight[j], tloop, right[j]);                 break;
	  case 2: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType2(loopRight[j], tloop, right[j], gamma1);         break;
	  case 3: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType3(loopRight[j], tloop, right[j]);                 break;
	  }
	}
	LOG(Message) << "loopRight packed for type " << type << std::endl;
	LOG(Message) << "Making EMF" << std::endl;
	A2ASpatialSum<SpinColourVector_v> spatial_sum;
	for ( unsigned int i = 0; i < N_i; i += cacheBlock )
	for ( unsigned int j = 0; j < N_j; j += cacheBlock ){
	  int Nii = MIN(N_i-i,cacheBlock);
	  int Njj = MIN(N_j-j,cacheBlock);

	  spatial_sum.Allocate(Nii, Njj, grid);
	  spatial_sum.PackLeft (leftv,     i, Nii);
	  spatial_sum.PackRight(loopRight, j, Njj);

	  Eigen::Tensor<ComplexD,3> emfBlock(nt, Nii, Njj);
	  emfBlock.setZero();
	  spatial_sum.Sum(emfBlock);

	  for(int t =0;t< nt;t++)
	  for(int ii=0;ii< Nii;ii++)
	  for(int jj=0;jj< Njj;jj++)
	    emf(0,0,t,i+ii,j+jj) = emfBlock(t,ii,jj);
	}// i,j
	LOG(Message) << "Done" << std::endl;

        double       blockSize, ioTime;

	std::string ioname = "type" + std::to_string(type) + "_" + nameg1_[ig] + "_" + nameg2_[ig];
	std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + ioname + ".h5";
        LOG(Message) << "Writing block to " << filename << std::endl;
        makeFileDir(filename, grid);

        //ioTime = -GET_TIMER("IO: write block");
        //START_TIMER("IO: total");
#ifdef HADRONS_A2AM_PARALLEL_IO
        grid->Barrier();
        LOG(Message) << "HADRONS_A2AM_PARALLEL_IO" << std::endl;
	if ( grid->ThisRank() == 0 ) {
#endif
	  // make task list for current node
	  A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename,ioname,nt, N_i, N_j);
	  A2AExtendedMesonFieldMetadata md;
	  md.gamma1 = nameg1_[ig];
	  md.gamma2 = nameg2_[ig];

	  // memory consuming
	  io.initFile(md, MAX(N_i,N_j));
	  //START_TIMER("IO: write block");
	  io.saveBlock(emf, 0, 0, 0, 0);
	  //STOP_TIMER("IO: write block");

#ifdef HADRONS_A2AM_PARALLEL_IO
	}
	grid->Barrier();
#endif
        //STOP_TIMER("IO: total");
      }// ig
    }// type

}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AExtendedMesonField_hpp_
