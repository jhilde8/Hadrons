/*
 * Test_writer_spread.cpp
 *
 * Reports what the per-timeslice file ownership rule actually produces on
 * this machine: whether every file has exactly one owner, and how far the
 * writers spread over ranks and nodes.
 *
 * Why this is a test and not instrumentation in the modules: the answer is a
 * pure function of the MPI decomposition, the node count, nt and the file
 * count. Nothing about the data, the trajectory or the hit count enters --
 * nFiles is nmom*ngamma (or ntypes*ngammaPairs) and ntOut is nt/Pt, both
 * fixed per job shape -- and OptimalCommunicator's relabeling is
 * deterministic for a given decomposition. So it needs running once per
 * (machine, --mpi, job shape), not once per job.
 *
 * Running it matters because index-space spread is NOT node-space spread:
 * ranks and coordinates are not related lexicographically once
 * OptimalCommunicator relabels for shared-memory locality (see
 * Grid/communicator/RingAllReduce.h, which records a wrong inverse from
 * assuming they were). The ownership arithmetic guarantees distinct ranks;
 * only measurement tells you about nodes.
 *
 * The rule mirrored here is the one in A2AMesonField, A2AExtendedMesonField
 * and A2AChromoMagneticOperatorField. Keep them in step.
 *
 * Usage:
 *   Test_writer_spread --grid 16.16.16.32 --mpi 1.1.2.4 [--ngroups 20]
 *
 * --ngroups is the number of file groups: nmom*ngamma for a meson field
 * (27 for the pion, 7 for the kaon), ntypes*ngammaPairs for an EMF (20),
 * nparities*nifOrthogs for a CMOF (4).
 */

#include <Grid/Grid.h>

using namespace Grid;

int main(int argc, char *argv[])
{
    Grid_init(&argc, &argv);

    unsigned int nFiles = 20;

    if (GridCmdOptionExists(argv, argv + argc, "--ngroups"))
    {
        nFiles = std::stoi(GridCmdOptionPayload(argv, argv + argc, "--ngroups"));
    }

    GridCartesian *grid = SpaceTimeGrid::makeFourDimGrid(
        GridDefaultLatt(), GridDefaultSimd(Nd, vComplex::Nsimd()),
        GridDefaultMpi());

    const int nd     = grid->Nd();
    const int nt     = grid->FullDimensions()[nd - 1];
    const int ntOut  = grid->LocalDimensions()[nd - 1];
    const int ct     = grid->ThisProcessorCoor()[nd - 1];
    const int Pt     = grid->ProcessorGrid()[nd - 1];
    const int nNodes = grid->NodeCount();
    const int nRank  = grid->RankCount();

    GRID_ASSERT(grid->LocalStarts()[nd - 1] == ct * ntOut);

    // Seat among the P_xyz ranks sharing this t coordinate, dimension 0
    // fastest; one loop yields both the seat and P_xyz.
    unsigned int mySeat = 0, P_xyz = 1;

    for (int mu = 0; mu < nd - 1; mu++)
    {
        mySeat += (unsigned int)grid->ThisProcessorCoor()[mu] * P_xyz;
        P_xyz  *= (unsigned int)grid->ProcessorGrid()[mu];
    }

    const unsigned int step = (P_xyz > (unsigned int)ntOut)
                            ? P_xyz / (unsigned int)ntOut : 1u;

    auto owns = [&](const unsigned int f, const int gt)
    {
        if (gt / ntOut != ct)
        {
            return false;
        }

        unsigned int base = (unsigned int)(((uint64_t)f * P_xyz) / nFiles);

        return (base + (unsigned int)(gt % ntOut) * step) % P_xyz == mySeat;
    };

    std::cout << GridLogMessage << "grid " << grid->FullDimensions()
              << "  mpi " << grid->ProcessorGrid()
              << "  shm " << grid->ShmGrid() << std::endl;
    std::cout << GridLogMessage << nRank << " ranks on " << nNodes
              << " nodes, P_xyz = " << P_xyz << ", P_t = " << Pt
              << ", nt = " << nt << ", ntOut = " << ntOut << std::endl;
    std::cout << GridLogMessage << nFiles << " file groups x " << nt
              << " timeslices = " << nFiles * nt << " files, seat step = "
              << step << std::endl;
    if (grid->ShmGrid()[nd - 1] != 1)
    {
        std::cout << GridLogMessage << "note: shm grid spans the time axis, so "
                  << "a node straddles " << grid->ShmGrid()[nd - 1]
                  << " t slabs and writers from different slabs can share a "
                  << "node. This costs node spread, not throughput -- the "
                  << "checksum is CPU work and co-resident ranks use "
                  << "different cores." << std::endl;
    }

    // owners[f*nt + gt] must sum to exactly 1 across the communicator.
    // nodeMark[f*nNodes + node] counts this node's writers of group f.
    // rankMark[f] counts distinct ranks writing group f.
    std::vector<uint64_t> owners((size_t)nFiles * nt, 0);
    std::vector<uint64_t> nodeMark((size_t)nFiles * nNodes, 0);
    std::vector<uint64_t> rankMark(nFiles, 0);
    const size_t          myNode = (size_t)GlobalSharedMemory::WorldNode;

    for (unsigned int f = 0; f < nFiles; ++f)
    {
        bool any = false;

        for (int gt = 0; gt < nt; ++gt)
        {
            if (owns(f, gt))
            {
                owners[(size_t)f * nt + gt] = 1;
                any = true;
            }
        }
        if (any)
        {
            nodeMark[(size_t)f * nNodes + myNode] = 1;
            rankMark[f]                           = 1;
        }
    }

    grid->GlobalSumVector(owners.data(), (int)owners.size());
    grid->GlobalSumVector(nodeMark.data(), (int)nodeMark.size());
    grid->GlobalSumVector(rankMark.data(), (int)nFiles);

    unsigned int nUnowned = 0, nMulti = 0;

    for (size_t k = 0; k < owners.size(); ++k)
    {
        if (owners[k] == 0) nUnowned++;
        if (owners[k] >  1) nMulti++;
    }

    unsigned int      minN = (unsigned int)nNodes, maxN = 0, minR = ~0u, maxR = 0;
    std::vector<char> nodeUsed(nNodes, 0);

    for (unsigned int f = 0; f < nFiles; ++f)
    {
        unsigned int nodes = 0;

        for (int n = 0; n < nNodes; ++n)
        {
            if (nodeMark[(size_t)f * nNodes + n])
            {
                nodes++;
                nodeUsed[n] = 1;
            }
        }
        minN = std::min(minN, nodes);
        maxN = std::max(maxN, nodes);
        minR = std::min(minR, (unsigned int)rankMark[f]);
        maxR = std::max(maxR, (unsigned int)rankMark[f]);
    }

    unsigned int nodesAll = 0;

    for (int n = 0; n < nNodes; ++n)
    {
        if (nodeUsed[n]) nodesAll++;
    }

    // Per group: instantaneous concurrency, and the NVMe pressure while one
    // group is being written. Across all groups: where the output ends up
    // sitting, which is what decides how parallel the stage-out can be.
    std::cout << GridLogMessage << "ranks/group " << minR << "-" << maxR
              << " (want " << std::min<unsigned int>(nt, nRank) << ")"
              << ", nodes/group " << minN << "-" << maxN
              << " (of " << nNodes << ")" << std::endl;
    std::cout << GridLogMessage << "nodes holding output across all groups: "
              << nodesAll << " of " << nNodes << std::endl;

    const bool ok = (nUnowned == 0) && (nMulti == 0);

    if (nUnowned)
    {
        std::cout << GridLogError << nUnowned << " file(s) have NO owner and "
                  << "would never be written" << std::endl;
    }
    if (nMulti)
    {
        std::cout << GridLogError << nMulti << " file(s) have MORE THAN ONE "
                  << "owner - concurrent writers would race" << std::endl;
    }
    std::cout << GridLogMessage << "ownership: " << (ok ? "OK" : "BROKEN")
              << std::endl;

    Grid_finalize();

    return ok ? 0 : 1;
}
