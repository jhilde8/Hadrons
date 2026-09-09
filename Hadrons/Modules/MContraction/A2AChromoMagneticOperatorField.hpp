/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid

Source file: Hadrons/Modules/MContraction/A2AChromoMagneticOperatorField.hpp

Copyright (C) 2015-2019

Author: Peter Boyle <paboyle@bnl.gov>
Author: Jonas Hildebrand <jonas.hildebrand@uconn.edu>
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
#include <Grid/algorithms/blas/A2ASpatialSum.h>
#include <iomanip>

BEGIN_HADRONS_NAMESPACE
/******************************************************************************
 *          All-to-all chromo-magnetic operator field creation (GPU)          *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class A2AChromoMagneticOperatorFieldPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(A2AChromoMagneticOperatorFieldPar,
                                    int,         block,
                                    int,         cacheBlock,
                                    std::string, parities,
                                    std::string, left,
                                    std::string, right,
                                    std::string, gauge,
                                    std::string, output,
                                    std::string, ifOrthogs,
                                    bool,        timeSliceIO);
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
  // One per module, not one per (ifOrthog, parity) pair -- see the same member
  // in TA2AExtendedMesonField. Nothing carries across a pair; a fresh object
  // cost 4 rounds of freeing and reallocating the largest device buffers here.
  A2ASpatialSum<vobj> spatial_sum_;
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
  auto &left    = envGet(std::vector<FermionField>, par().left);
  auto &right   = envGet(std::vector<FermionField>, par().right);
  const auto &U = envGet(GaugeField, par().gauge);

  GridBase *grid = left[0].Grid();

  LOG(Message) << "Computing all-to-all ChromoMagnetic operator fields (GPU)" << std::endl;

  int nt         = env().getDim().back();
  int N_i        = left.size();
  int N_j        = right.size();
  int block = par().block;
  int cacheBlock = par().cacheBlock;

  // timeSliceIO: SumRing stops after the spatial reduce, so this rank holds
  // only its own t slab and writes one file per timeslice it owns. Global
  // timeslice of local index lt is ct*ntOut + lt.
  //
  // There is no staging buffer: each (i,j) block is written as a hyperslab
  // straight out of SumRing's result tensor, the way A2AMesonField and
  // A2AFewMesonField already do it. Buffering the whole field first cost
  // ntOut*N_i*N_j per rank, which scales with the square of the hit count.
  const bool tsIO = par().timeSliceIO;
  int nd     = grid->Nd();
  int ct     = grid->ThisProcessorCoor()[nd - 1];
  int ntOut  = tsIO ? grid->LocalDimensions()[nd - 1] : nt;
  int ntFile = tsIO ? 1 : nt;

  // Seat among the P_xyz ranks sharing this t coordinate. The (ifOrthog,
  // parity) loop is sequential, so all concurrency comes from the timeslice
  // axis; see the header comment.
  unsigned int mySeat = 0, P_xyz = 1;
  for (int mu = 0; mu < nd - 1; mu++)
  {
      mySeat += (unsigned int)grid->ThisProcessorCoor()[mu] * P_xyz;
      P_xyz  *= (unsigned int)grid->ProcessorGrid()[mu];
  }

  unsigned int nFiles = (unsigned int)(ifOrthogs_.size() * parities_.size());
  auto ownerFn = [tsIO, ntOut, ct, P_xyz, mySeat, nFiles, grid]
                 (const unsigned int f, const int gt)
  {
      if (!tsIO)
          return grid->ThisRank() == 0;
      if (gt / ntOut != ct)
          return false;

      unsigned int base = (unsigned int)(((uint64_t)f * P_xyz) / nFiles);
      unsigned int step = (P_xyz > (unsigned int)ntOut)
                        ? P_xyz / (unsigned int)ntOut : 1u;
      return (base + (unsigned int)(gt % ntOut) * step) % P_xyz == mySeat;
  };
  unsigned int fileIdx = 0;

  LOG(Message) << "Left: '"        << par().left  << "' Right: '"
               << par().right      << "'"          << std::endl;
  LOG(Message) << "Gauge field: '" << par().gauge  << "'"          << std::endl;
  LOG(Message) << "Parities:"      << std::endl;
  for (auto &p: parities_)
    LOG(Message) << "  " << p << std::endl;
  LOG(Message) << "CMO field size: " << nt << "*" << N_i << "*" << N_j
               << " (filesize " << sizeString(ntFile*N_i*N_j*sizeof(HADRONS_A2AM_IO_TYPE))
               << (tsIO ? "/timeslice)" : ")") << std::endl;
  if (tsIO)
      LOG(Message) << "Per-timeslice IO: this rank holds t = " << ct*ntOut
                   << ".." << ct*ntOut + ntOut - 1 << std::endl;

  std::vector<FermionField> loopRight(block, grid);

  std::array<double, 6> sumTimings = {};
  std::array<double, 6> sumBytes   = {};
  std::array<double, 7> ioTimings  = {};

  // Every rank makes the output directory itself rather than relying on
  // makeFileDir's boss-only mkdir: under timeSliceIO the writers are spread
  // across ranks, and on node-local storage (an NVMe burst buffer) a directory
  // made on the boss node is simply absent on every other one. Checked,
  // because a discarded failure resurfaces as an opaque "errno = 2" out of
  // H5Fcreate much further downstream.
  std::string dirBase = par().output + "." + std::to_string(vm().getTrajectory());

  startTimer("mkdir");
  if (Hadrons::mkdir(dirBase))
  {
      HADRONS_ERROR(Io, "cannot create directory '" + dirBase + "' ("
                        + std::strerror(errno) + ")");
  }
  grid->Barrier();
  stopTimer("mkdir");

  // Regenerated per write rather than cached: A2AMatrixIo holds only names and
  // dimensions, never an open handle -- saveBlock opens and closes itself -- so
  // rebuilding one costs two string copies against a 256 KiB write.
  auto filenameFn = [&dirBase, tsIO](const std::string &ioname, const int gt)
  {
      std::stringstream fn;

      fn << dirBase << "/" << ioname;
      if (tsIO) fn << ".t" << std::setfill('0') << std::setw(4) << gt;
      fn << ".h5";

      return fn.str();
  };

  for (auto &ifOrthog: ifOrthogs_) {
    std::vector<GaugeMat>  G;
    Vector<Gamma::Algebra> Sigma;

    startTimer("CMO contraction");
    if (ifOrthog == 0)
      Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContraction0(G, Sigma, U);
    else
      Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContraction1(G, Sigma, U);
    stopTimer("CMO contraction");
    LOG(Message) << "Field strength constructed for ifOrthog=" << ifOrthog << std::endl;

    for (auto &parity: parities_) {
      LOG(Message) << "Starting calculation with ifOrthog=" << ifOrthog
                   << " parity=" << parity << std::endl;
      std::string ioname = "parity" + std::to_string(parity);
      if (ifOrthog == 1)
        ioname = ioname + "_GijSij";
      else
        ioname = ioname + "_GitSit";

      // Create the timeslice files this rank owns before the block sweep, then
      // write each (i,j) block into them as it is produced. Both loops test the
      // same ownerFn(fileIdx, gt), so the rank that creates a file is always
      // the rank that writes it -- creating on one rank and writing from
      // another races on client-side metadata caching, since a plain Barrier
      // does not make the file visible from another node. The chunk is `block`,
      // matching the write granularity, so a hyperslab covers exactly one chunk
      // and never forces a read-modify-write.
      A2AChromoMagneticOperatorFieldMetadata md;
      md.meta = ioname;
      int nOwned = 0;
      for (int lt = 0; lt < (tsIO ? ntOut : 1); ++lt)
      {
        int gt = tsIO ? ct*ntOut + lt : 0;

        if (!ownerFn(fileIdx, gt)) continue;

        A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filenameFn(ioname, gt), ioname,
                                             ntFile, N_i, N_j);
        startTimer("initFile");
        io.initFile(md, block);
        stopTimer("initFile");
        nOwned++;
      }

      LOG(Message) << "Writing " << (tsIO ? nt : 1) << " file(s) to "
                   << dirBase << "/" << ioname << std::endl;
#ifdef HADRONS_A2AM_PARALLEL_IO
      startTimer("Barrier");
      grid->Barrier();
      stopTimer("Barrier");
#endif
      double writeTime = 0.;

      // Result buffers, one per distinct block shape. A block is full or on
      // the tail in each axis independently, so a 2x2 pool indexed by (i on
      // tail, j on tail) covers every case. Each slot is asked for the same
      // dimensions every time it is selected, so the first visit allocates and
      // every later one is a dimension assignment -- Eigen's resize
      // reallocates only when the total element count changes. That replaces
      // one construct/destruct per (i,j) block of a buffer SumRing overwrites
      // in full anyway.
      //
      // RowMajor is what lets SumRing take its direct device->host path: the
      // gathered panel's [gt][iii][m][jjj] layout and a RowMajor (ntOut, Nii, 1,
      // Njj) tensor are then the same addresses, so its scatter is skipped and
      // its "scatter" timer stays at zero. ColMajor would put t fastest in
      // memory while the copy-out below walks t outermost, which is both the
      // wrong order for that loop and the reason the direct path could not
      // apply. Element access is layout independent, so the values are
      // unchanged.
      Eigen::Tensor<ComplexD, 4, Eigen::RowMajor> resPool[2][2];

      LOG(Message) << "Making CMF" << std::endl;

      for (unsigned int j = 0; j < N_j; j += block) {
        int Njj = MIN(N_j-j, block);

        startTimer("CMOContractRight");
        for (int jj = 0; jj < Njj; jj++)
          Grid::A2AChromoMagneticOperator<GImpl,FImpl>::CMOContractRight(
              loopRight[jj], G, Sigma, right[j+jj], parity);
        stopTimer("CMOContractRight");

        startTimer("Allocate");
        spatial_sum_.AllocateRight(Njj, grid);
        stopTimer("Allocate");
        startTimer("Pack vectors");
        spatial_sum_.PackRight(loopRight, 0, Njj);
        stopTimer("Pack vectors");

        for (unsigned int i = 0; i < N_i; i += block) {
          int Nii = MIN(N_i-i, block);

          startTimer("Allocate");
          spatial_sum_.AllocateLeft(Nii);
          stopTimer("Allocate");
          startTimer("Pack vectors");
          spatial_sum_.PackLeftConj(left, i, Nii);
          stopTimer("Pack vectors");

          // Rank 4 with a singleton momentum axis: SumRing writes
          // result[t][i][m][j] for the general nmom case, and the CMO field
          // carries no momentum projection.
          //
          // No setZero: SumRing writes every element of the tensor on both its
          // direct and its scatter path, so zeroing first is dead work.
          auto &cmfBlock = resPool[Nii != block][Njj != block];
          startTimer("Allocate");
          cmfBlock.resize(ntOut, Nii, 1, Njj);
          stopTimer("Allocate");

          startTimer("Sum");
          spatial_sum_.SumRing(cmfBlock, cacheBlock, &sumTimings, &sumBytes, tsIO);
          stopTimer("Sum");

          // Straight out of the result tensor: cmfBlock is RowMajor
          // (ntOut, Nii, 1, Njj), so the Nii x Njj slab at fixed t is
          // contiguous and needs no staging. Without tsIO there is one file
          // with ntFile = nt = ntOut, and saveBlock's count of {nt, Nii, Njj}
          // spans the whole tensor from lt = 0.
          startTimer("IO");
          double dt = -usecond();
          for (int lt = 0; lt < (tsIO ? ntOut : 1); ++lt)
          {
            int gt = tsIO ? ct*ntOut + lt : 0;

            if (!ownerFn(fileIdx, gt)) continue;

            A2AMatrixIo<HADRONS_A2AM_IO_TYPE> io(filenameFn(ioname, gt), ioname,
                                                 ntFile, N_i, N_j);
            io.saveBlock(&cmfBlock(lt, 0, 0, 0), i, j, Nii, Njj, "", &ioTimings);
          }
          dt += usecond();
          writeTime += dt;
          stopTimer("IO");

          //LOG(Message) << "CMF made for i-block " << i/block
          //             << " j-block "             << j/block
          //             << " ifOrthog="            << ifOrthog
          //             << " parity="              << parity << std::endl;

        }// i
      }// j

      LOG(Message) << "CMF made for ifOrthog=" << ifOrthog << " parity=" << parity << std::endl;

      fileIdx++;
#ifdef HADRONS_A2AM_PARALLEL_IO
      startTimer("Barrier");
      grid->Barrier();
      stopTimer("Barrier");
#endif
      // Per rank, not global: writeTime accumulates this rank's own saveBlock
      // calls and ioBytes counts only the files it owns, so the rate is what
      // this rank achieved rather than the whole field over a barrier window.
      double ioBytes = static_cast<double>(nOwned) * ntFile * N_i * N_j
                       * sizeof(HADRONS_A2AM_IO_TYPE);
      if (writeTime > 0.)
          LOG(Message) << "IO ifOrthog=" << ifOrthog << " parity=" << parity
                       << ": " << sizeString(ioBytes) << " in " << writeTime
                       << " us local (" << ioBytes / writeTime * 1.e6 / 1024. / 1024.
                       << " MB/s effective)" << std::endl;
    }// parity
  }// ifOrthog

  // Throughput of the post-GEMM SumRing stages -- bytesMoved[k] and
  // sumTimings[k] accumulate the same way across all (ifOrthog,parity,i,j)
  // calls and all cacheBlock tiles, so their ratio is the average effective
  // bandwidth of that stage over the whole run, comparable across different
  // cacheBlock choices.
  //
  // The two ring stages report bytes on the wire rather than payload, so
  // their rates are the ones comparable with a link rate; the local stages
  // report the bytes they actually touch. See the SumRing header comment.
  auto gbps = [](double bytes, double us)
  {
      return (us > 0.) ? bytes / us * 1.e6 / 1024. / 1024. / 1024. : 0.;
  };
  LOG(Message) << "Sum detail (us), rank 0:" << std::endl;
  LOG(Message) << "  GEMM            = " << sumTimings[0] << std::endl;
  LOG(Message) << "  device->host    = " << sumTimings[1]
               << " (" << gbps(sumBytes[1], sumTimings[1]) << " GB/s)" << std::endl;
  LOG(Message) << "  gather to slab  = " << sumTimings[2]
               << " (" << gbps(sumBytes[2], sumTimings[2]) << " GB/s)" << std::endl;
  LOG(Message) << "  spatial reduce  = " << sumTimings[3]
               << " (" << gbps(sumBytes[3], sumTimings[3]) << " GB/s wire)" << std::endl;
  LOG(Message) << "  scatter         = " << sumTimings[4]
               << " (" << gbps(sumBytes[4], sumTimings[4]) << " GB/s)" << std::endl;
  LOG(Message) << "  temporal gather = " << sumTimings[5]
               << " (" << gbps(sumBytes[5], sumTimings[5]) << " GB/s wire)" << std::endl;
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

#endif // Hadrons_MContraction_A2AChromoMagneticOperatorField_hpp_
