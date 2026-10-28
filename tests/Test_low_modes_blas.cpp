// Workaround for NVCC 11.5 + GCC 11 std::function incompatibility.
// cudafe mishandles _ArgTypes... in std::function's constrained constructor.
#ifdef __CUDACC__
#pragma nv_diag_suppress 20011
#endif

/*
 * Test_low_modes_blas.cpp
 *
 * In-memory regression of the GEMM promotion in
 * MUtilities::A2ALowModesCoarseBLAS against the blockPromote of
 * MUtilities::A2ALowModesCoarse, on the first --size coarse vectors of the
 * 64I light eigenpack (EPACK_STEM, 200 fine) of configuration --traj
 * (required). Nothing is written. The pack is
 * built on Hadrons' own environment grids, the ones
 * MIO::LoadCoarseFermionEigenPack200F uses, so the blocking is the production
 * one.
 *
 * Correctness. The modes are promoted in batches of --batch with
 * MultiRHSBlockProject, as the module does, and every mode is compared with
 * a scalar blockPromote of the same coarse vector. Pass is
 * |blas - scalar|/|scalar| <= --tol for every mode: a GEMM and a loop of axpy
 * sum in different orders, so the two are not bit-identical. Both stay in
 * single precision; the precisionChange after them is the same in both
 * modules and is not tested. Pick --size so the last batch is short, as in
 * production.
 *
 * Timing. Each --sweep entry nrhs promotes the first nrhs modes in one call,
 * after one untimed call that loads rocBLAS's kernels, and reports the time
 * per call and per mode next to the scalar time per mode. Entries larger than
 * --size are skipped.
 *
 * Usage:
 *   srun ... ./Test_low_modes_blas --grid 64.64.64.128 --mpi 4.4.4.8 \
 *       --traj 1540 --size 130 --batch 64 --sweep 1,8,16,32,64,128
 */

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>
#include "TestMFShells.hpp"

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

static void promoteBatch(MultiRHSBlockProject<FermionFieldF> &blockProj, EPack &epack,
                         const unsigned int first, const unsigned int n,
                         std::vector<FermionFieldF> &fine)
{
    std::vector<CoarseField> coarse(n, epack.evecCoarse[0].Grid());

    for (unsigned int j = 0; j < n; j++)
    {
        coarse[j] = epack.evecCoarse[first + j];
    }
    blockProj.blockPromote(fine, coarse);
    for (auto &f: fine)
    {
        f.Checkerboard() = epack.evec[0].Checkerboard();
    }
}

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
    int          traj      = 0;
    unsigned int size      = 130;
    unsigned int batch     = 64;
    double       tol       = 1.0e-5;
    std::string  sweepList = "1,8,16,32,64,128";

    if (GridCmdOptionExists(argv, argv + argc, "--traj"))
        traj      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--traj"));
    if (GridCmdOptionExists(argv, argv + argc, "--size"))
        size      = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--size"));
    if (GridCmdOptionExists(argv, argv + argc, "--batch"))
        batch     = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--batch"));
    if (GridCmdOptionExists(argv, argv + argc, "--tol"))
        tol       = std::stod(GridCmdOptionPayload(argv, argv + argc, "--tol"));
    if (GridCmdOptionExists(argv, argv + argc, "--sweep"))
        sweepList = GridCmdOptionPayload(argv, argv + argc, "--sweep");

    std::vector<unsigned int> sweep = strToVec<unsigned int>(cliListToPar(sweepList));

    if (!GridCmdOptionExists(argv, argv + argc, "--traj"))
    {
        HADRONS_ERROR(Argument, "--traj is required: the configuration of the eigenpack");
    }

    // ------------------------------------------------------------------
    // Eigenpack on the environment grids, and the packed basis.
    // ------------------------------------------------------------------
    auto     &env    = Environment::getInstance();
    GridBase *gridRb = env.getRbGrid<FermionFieldF::vector_type>(LS);
    GridBase *gridC  = env.getCoarseGrid<CoarseField::vector_type>(BLOCK_SIZE, LS);
    EPack    epack(200, size, gridRb, gridC);

    epack.read(EPACK_STEM, true, traj);

    MultiRHSBlockProject<FermionFieldF> blockProj;
    double                              t = -usecond();

    blockProj.Allocate(epack.evec.size(), gridRb, gridC);
    blockProj.ImportBasis(epack.evec);
    t += usecond();
    LOG(Message) << "Packed " << epack.evec.size() << " fine basis vectors in "
                 << t/1.0e6 << " s, " << blockProj.DeviceBytes()/1024./1024./1024.
                 << " GiB resident device memory per rank" << std::endl;

    // ------------------------------------------------------------------
    // Correctness.
    // ------------------------------------------------------------------
    FermionFieldF ref(gridRb), diff(gridRb);
    double        worst = 0., tScalar = 0.;
    unsigned int  worstMode = 0;

    ref.Checkerboard() = epack.evec[0].Checkerboard();
    for (unsigned int il0 = 0; il0 < size; il0 += batch)
    {
        unsigned int               n = std::min(batch, size - il0);
        std::vector<FermionFieldF> fine(n, gridRb);
        double                     batchWorst = 0.;

        promoteBatch(blockProj, epack, il0, n, fine);
        for (unsigned int j = 0; j < n; j++)
        {
            tScalar -= usecond();
            blockPromote(epack.evecCoarse[il0 + j], ref, epack.evec);
            tScalar += usecond();

            diff = ref - fine[j];

            double rel = std::sqrt(norm2(diff)/norm2(ref));

            batchWorst = std::max(batchWorst, rel);
            if (rel > worst)
            {
                worst     = rel;
                worstMode = il0 + j;
            }
        }
        LOG(Message) << "Modes " << il0 << "-" << il0 + n - 1
                     << ": max relative difference " << batchWorst << std::endl;
    }
    tScalar /= size;

    // ------------------------------------------------------------------
    // Timing sweep.
    // ------------------------------------------------------------------
    for (auto nrhs: sweep)
    {
        if (nrhs > size)
        {
            LOG(Warning) << "Sweep entry " << nrhs << " > --size " << size
                         << ", skipped" << std::endl;
            continue;
        }

        std::vector<FermionFieldF> fine(nrhs, gridRb);

        promoteBatch(blockProj, epack, 0, nrhs, fine);
        t = -usecond();
        promoteBatch(blockProj, epack, 0, nrhs, fine);
        t += usecond();
        LOG(Message) << "nrhs " << nrhs << ": " << t/1.0e3 << " ms per call, "
                     << t/1.0e3/nrhs << " ms per mode (scalar " << tScalar/1.0e3
                     << " ms per mode)" << std::endl;
    }

    bool pass = (worst <= tol);

    LOG(Message) << "Max relative difference " << worst << " at mode " << worstMode
                 << ", tolerance " << tol << ": " << (pass ? "PASS" : "FAIL")
                 << std::endl;

    Grid_finalize();

    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
