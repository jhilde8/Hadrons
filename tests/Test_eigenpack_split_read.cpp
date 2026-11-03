// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_eigenpack_split_read.cpp
 *
 * Regression of the split-grid coarse read of MIO::LoadCoarseEigenPackSplit
 * (CoarseEigenPack::readCoarseSplit) against the full-grid readCoarse, on the
 * first --size coarse vectors of the 64I light eigenpack (EPACK_STEM) of
 * configuration --traj (required). The packs are built on Hadrons' own
 * environment grids, as MIO::LoadCoarseFermionEigenPack200F does, and the
 * split grid is built as in the module from the 4D layout --io-split.
 *
 * Pass is bitwise: the read is pure data movement, so every coarse vector
 * must satisfy norm2(split - full) == 0 and every eval must be equal.
 *
 * Timing. Both reads are timed, and the split one logs its per-partition
 * timings. The split read runs first so it is never helped by a warm client
 * cache; the full-grid time that follows is, if anything, flattered. Pick
 * --size so the last batch is short, as in production.
 *
 * Usage:
 *   srun ... ./Test_eigenpack_split_read --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --traj 1540 --io-split 1.1.2.4 --size 2000
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

typedef FIMPLF::FermionField                   FermionFieldF;
typedef CoarseFermionEigenPack<FIMPLF, 200>    EPack;
typedef EPack::CoarseField                     CoarseField;

// PyXmlGen/config.py CG_EIGENPACK, as in Test_split_write.
static const std::string EPACK_STEM =
    "/lustre/orion/phy157/proj-shared/phy157_dwf/jhilde/main_64I/converted/vec";
static const unsigned int     LS         = 12;
static const std::vector<int> BLOCK_SIZE = {4, 4, 4, 4, 12};

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());

    // ------------------------------------------------------------------
    // Parse CLI arguments.
    // ------------------------------------------------------------------
    int          traj    = 0;
    unsigned int size    = 2000;
    std::string  ioSplit = "1.1.2.4";

    if (!GridCmdOptionExists(argv, argv + argc, "--traj"))
    {
        HADRONS_ERROR(Argument, "--traj is required: the configuration of the eigenpack");
    }
    traj = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--traj"));
    if (GridCmdOptionExists(argv, argv + argc, "--size"))
        size    = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--size"));
    if (GridCmdOptionExists(argv, argv + argc, "--io-split"))
        ioSplit = GridCmdOptionPayload(argv, argv + argc, "--io-split");

    // ------------------------------------------------------------------
    // Packs on the environment grids, and the split coarse grid.
    // ------------------------------------------------------------------
    auto          &env   = Environment::getInstance();
    GridBase      *gridRb = env.getRbGrid<FermionFieldF::vector_type>(LS);
    GridCartesian *cg     = env.getCoarseGrid<CoarseField::vector_type>(BLOCK_SIZE, LS);
    GridCartesian *cg4d   = env.getCoarseGrid<CoarseField::vector_type>(BLOCK_SIZE);
    EPack         packSplit(0, size, gridRb, cg), packFull(0, size, gridRb, cg);
    Coordinate    splitMpi;

    GridCmdOptionIntVector(ioSplit, splitMpi);

    GridCartesian                  split4d(cg4d->FullDimensions(), cg4d->_simd_layout,
                                           splitMpi, *cg4d);
    std::unique_ptr<GridCartesian> split(SpaceTimeGrid::makeFiveDimGrid(
                                             cg->FullDimensions()[0], &split4d));
    unsigned int                   nPart = cg->ProcessorCount()/split->ProcessorCount();
    double                         bytes = static_cast<double>(cg->gSites())
                                           *sizeof(CoarseField::scalar_object);
    std::vector<double>            t;
    double                         tSplit, tFull;

    LOG(Message) << "Split coarse read: " << nPart << " partitions of layout "
                 << ioSplit << std::endl;

    // ------------------------------------------------------------------
    // Reads, split first.
    // ------------------------------------------------------------------
    tSplit = -usecond();
    packSplit.readCoarseSplit(EPACK_STEM, split.get(), t, traj);
    tSplit += usecond();
    A2AVectorsIo::logPartitions("Read coarse", 0, size, nPart, t, bytes);

    tFull = -usecond();
    packFull.readCoarse(EPACK_STEM, true, traj);
    tFull += usecond();

    // ------------------------------------------------------------------
    // Bitwise comparison.
    // ------------------------------------------------------------------
    CoarseField  diff(cg);
    unsigned int nBad = 0;

    for (unsigned int k = 0; k < size; ++k)
    {
        diff = packSplit.evecCoarse[k] - packFull.evecCoarse[k];

        double d = norm2(diff);

        if (d != 0. || packSplit.evalCoarse[k] != packFull.evalCoarse[k])
        {
            LOG(Message) << "Vector " << k << ": norm2(diff) " << d << ", eval split "
                         << packSplit.evalCoarse[k] << " full "
                         << packFull.evalCoarse[k] << std::endl;
            nBad++;
        }
    }

    bool pass = (nBad == 0);

    LOG(Message) << "Split read " << tSplit/1.0e6 << " s, full-grid read "
                 << tFull/1.0e6 << " s, speedup " << tFull/tSplit << std::endl;
    LOG(Message) << nBad << " of " << size << " coarse vectors differ: "
                 << (pass ? "PASS" : "FAIL") << std::endl;

    Grid_finalize();

    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
