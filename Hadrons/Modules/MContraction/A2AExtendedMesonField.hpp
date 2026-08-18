/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Grid/algorithms/multigrid/GeneralCoarsenedMatrix.hGrid/algorithms/multigrid/GeneralCoarsenedMatrix.hSource file: Hadrons/Modules/MContraction/A2AExtendedMesonField.hpp

Copyright (C) 2015-2019

Author: Peter Boyle <paboyle@bnl.gov>
Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
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
                                    int, block,
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
    envCreateLat(PropagatorField, getName() + "_propagatorLoop");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TA2AExtendedMesonField<FImpl>::execute(void)
{
  typedef iSpinColourVector<vector_type> SpinColourVector_v;
    auto &left  = envGet(std::vector<FermionField>, par().left);
    auto &right = envGet(std::vector<FermionField>, par().right);
    auto &loop1 = envGet(std::vector<const FermionField *>, par().loop_vw1);
    auto &loop2 = envGet(std::vector<const FermionField *>, par().loop_vw2);

    GridBase *grid = left[0].Grid();

    LOG(Message) << "Computing all-to-all EXTENDED meson fields" << std::endl;

    int nt         = env().getDim().back();
    int N_i        = left.size();
    int N_j        = right.size();
    int block = par().block;
    int cacheBlock = par().cacheBlock;
    Vector<HADRONS_A2AM_IO_TYPE> mBuf; mBuf.resize(nt*N_i*N_j);

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

    auto &loop = envGet(PropagatorField, getName() + "_propagatorLoop");
    startTimer("LoopPropagator");
    Grid::A2AExtendedMesonField<FImpl>::LoopPropagatorPtr(loop, loop1, loop2);
    stopTimer("LoopPropagator");
    LOG(Message) << "Quark loop calculated" << std::endl;

    std::vector<FermionField> loopRight(block, grid);
    PropagatorField tloop(grid);

    std::array<double, 5> sumTimings = {};
    std::array<double, 5> sumBytes   = {};
    std::array<double, 7> ioTimings  = {};

    for (int &type: types_){

      for (int ig = 0 ; ig < gamma1_.size() ; ++ig ){

	A2AMatrixSet<HADRONS_A2AM_IO_TYPE> emf(mBuf.data(),1,1,nt,N_i,N_j);

	Vector<Gamma::Algebra> gamma1(gamma1_[ig].begin(), gamma1_[ig].end());
	Vector<Gamma::Algebra> gamma2(gamma2_[ig].begin(), gamma2_[ig].end());

	startTimer("Loop contraction");
	tloop = Zero();
	switch (type) {
	case 0: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType0(tloop, loop);                 break;
	case 1: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType1(tloop, loop, gamma1, gamma2); break;
	case 2: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType2(tloop, loop, gamma2);         break;
	case 3: Grid::A2AExtendedMesonField<FImpl>::LoopContractionType3(tloop, loop, gamma1, gamma2); break;
	}
	stopTimer("Loop contraction");
	//LOG(Message) << "tloop contraction done for type " << type << std::endl;
	//LOG(Message) << "Making EMF" << std::endl;
	A2ASpatialSum<SpinColourVector_v> spatial_sum;

	for ( unsigned int j = 0; j < N_j; j += block ){
	  int Njj = MIN(N_j-j,block);
	  startTimer("LoopRight contraction");
	  for (int jj = 0; jj < Njj; jj++) {
	    switch (type) {
	    case 0: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType0(loopRight[jj], tloop, right[j+jj], gamma1, gamma2); break;
	    case 1: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType1(loopRight[jj], tloop, right[j+jj]);                 break;
	    case 2: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType2(loopRight[jj], tloop, right[j+jj], gamma1);         break;
	    case 3: Grid::A2AExtendedMesonField<FImpl>::LoopRightContractionType3(loopRight[jj], tloop, right[j+jj]);                 break;
	    }
	  }
	  stopTimer("LoopRight contraction");
	  //LOG(Message) << "loopRight packed for j-block " << j/block << " type " << type << std::endl;

	  startTimer("Allocate");
	  spatial_sum.AllocateRight(Njj, grid);
	  stopTimer("Allocate");
	  startTimer("Pack vectors");
	  spatial_sum.PackRight(loopRight, 0, Njj);
	  stopTimer("Pack vectors");

	  for ( unsigned int i = 0; i < N_i; i += block ) {
	    int Nii = MIN(N_i-i,block);

	    startTimer("Allocate");
	    spatial_sum.AllocateLeft(Nii);
	    stopTimer("Allocate");
	    startTimer("Pack vectors");
	    spatial_sum.PackLeftConj(left, i, Nii);
	    stopTimer("Pack vectors");

	    Eigen::Tensor<ComplexD,3> emfBlock(nt, Nii, Njj);
	    emfBlock.setZero();
	    startTimer("Sum");
	    // spatial_sum.Sum(emfBlock, &sumTimings, &sumBytes);
	    spatial_sum.SumCacheBlocked(emfBlock, cacheBlock, &sumTimings, &sumBytes);
	    stopTimer("Sum");

	    for(int t =0;t< nt;t++)
	      for(int ii=0;ii< Nii;ii++)
		for(int jj=0;jj< Njj;jj++)
		  emf(0,0,t,i+ii,j+jj) = emfBlock(t,ii,jj);
	}
	//LOG(Message) << "EMF made for j-block " << j/block << " type " << type << std::endl;

	}// i,j
	LOG(Message) << "EMF made for type " << type << "; gamma1: " << nameg1_[ig] << "; gamma2: " << nameg2_[ig] << std::endl;

	std::string ioname = "type" + std::to_string(type) + "_" + nameg1_[ig] + "_" + nameg2_[ig];
	std::string filename = par().output + "." + std::to_string(vm().getTrajectory()) + "/" + ioname + ".h5";
        LOG(Message) << "Writing block to " << filename << std::endl;
        makeFileDir(filename, grid);
        double ioBytes = static_cast<double>(nt) * N_i * N_j * sizeof(HADRONS_A2AM_IO_TYPE);
        startTimer("IO");
        double writeTime = -usecond();
#ifdef HADRONS_A2AM_PARALLEL_IO
        startTimer("Barrier");
        grid->Barrier();
        stopTimer("Barrier");
        LOG(Message) << "HADRONS_A2AM_PARALLEL_IO" << std::endl;
	if ( grid->ThisRank() == 0 ) {
#endif
	  A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filename,ioname,nt, N_i, N_j);
	  A2AExtendedMesonFieldMetadata md;
	  md.gamma1 = nameg1_[ig];
	  md.gamma2 = nameg2_[ig];
	  startTimer("initFile");
	  io.initFile(md, MAX(N_i,N_j));
	  stopTimer("initFile");
	  io.saveBlock(emf, 0, 0, 0, 0, &ioTimings);
#ifdef HADRONS_A2AM_PARALLEL_IO
	}
	startTimer("Barrier");
	grid->Barrier();
	stopTimer("Barrier");
#endif
        writeTime += usecond();
        stopTimer("IO");
        // writeTime is this rank's own wall time from barrier to barrier
        // (only rank 0's LOG output survives, since Grid_quiesce_nodes
        // suppresses the rest by default).
        if (writeTime > 0.)
            LOG(Message) << "IO type=" << type << " ig=" << ig << ": "
                         << sizeString(ioBytes) << " in " << writeTime
                         << " us local (" << ioBytes / writeTime * 1.e6 / 1024. / 1024.
                         << " MB/s effective)" << std::endl;
      }// ig
    }// type

    // Throughput of the host-side, post-GEMM Sum() stages -- bytesMoved[k]
    // and sumTimings[k] accumulate the same way across all (type,ig,i,j)
    // calls and all cacheBlock tiles, so their ratio is the average
    // effective bandwidth of that stage over the whole run, comparable
    // across different cacheBlock choices.
    auto gbps = [](double bytes, double us)
    {
        return (us > 0.) ? bytes / us * 1.e6 / 1024. / 1024. / 1024. : 0.;
    };
    LOG(Message) << "Sum detail (us), rank 0:" << std::endl;
    LOG(Message) << "  GEMM            = " << sumTimings[0] << std::endl;
    LOG(Message) << "  device->host    = " << sumTimings[1]
                 << " (" << gbps(sumBytes[1], sumTimings[1]) << " GB/s)" << std::endl;
    LOG(Message) << "  transpose-1     = " << sumTimings[2]
                 << " (" << gbps(sumBytes[2], sumTimings[2]) << " GB/s)" << std::endl;
    LOG(Message) << "  GlobalSumVector = " << sumTimings[3]
                 << " (" << gbps(sumBytes[3], sumTimings[3]) << " GB/s)" << std::endl;
    LOG(Message) << "  transpose-2     = " << sumTimings[4]
                 << " (" << gbps(sumBytes[4], sumTimings[4]) << " GB/s)" << std::endl;
    LOG(Message) << "IO detail (us), rank 0:" << std::endl;
    LOG(Message) << "  open            = " << ioTimings[0]  << std::endl;
    LOG(Message) << "  push/group      = " << ioTimings[1]  << std::endl;
    LOG(Message) << "  openDataSet     = " << ioTimings[2]  << std::endl;
    LOG(Message) << "  getSpace        = " << ioTimings[3]  << std::endl;
    LOG(Message) << "  selectHyperslab = " << ioTimings[4]  << std::endl;
    LOG(Message) << "  write           = " << ioTimings[5]  << std::endl;
    LOG(Message) << "  close(fsync)    = " << ioTimings[6]  << std::endl;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_A2AExtendedMesonField_hpp_
