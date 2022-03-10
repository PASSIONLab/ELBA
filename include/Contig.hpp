
#ifndef CONTIG_HPP
#define CONTIG_HPP

#include <cmath>
#include <map>
#include <fstream>

#include <cassert>

#include "TraceUtils.hpp"
#include "ReadOverlap.hpp"
#include "Utils.hpp"
#include "CC.h"

using namespace combblas;

FullyDistVec<int64_t, int64_t> GetContigAssignments
(
    const SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>& OverlapGraph,
    FullyDistVec<int64_t, int64_t>& Branches,
    int64_t& NumContigs
)
{
    SpParMat<int64_t, int64_t, SpDCCols<int64_t, int64_t>> A = OverlapGraph;

    SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> D1 = A;
    FullyDistVec<int64_t, int64_t> degs1(D1.getcommgrid());

    D1.Reduce(degs1, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    Branches = degs1.FindInds(bind2nd(std::greater<int64_t>(), 2));

    A.PruneFull(Branches, Branches);

    SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> D2 = A;
    FullyDistVec<int64_t, int64_t> degs2(D2.getcommgrid());

    D2.Reduce(degs2, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    FullyDistVec<int64_t, int64_t> vCC = CC(A, NumContigs);

    return vCC;
}

FullyDistVec<int64_t, int64_t> GetContigSizes
(
    const FullyDistVec<int64_t, int64_t>& Assignments,
    const int64_t NumContigs
)
{
    std::vector<int64_t> LocalCCSizes(NumContigs, 0);
    std::vector<int64_t> LocalCC = Assignments.GetLocVec();

    int64_t LocalCCSize = Assignments.LocArrSize();

    for (int64_t i = 0; i < LocalCCSize; ++i)
        LocalCCSizes[LocalCC[i]]++;

    int nprocs, myrank;

    MPI_Comm World = Assignments.getcommgrid()->GetWorld();
    MPI_Comm_size(World, &nprocs);
    MPI_Comm_rank(World, &myrank);

    int avesize = NumContigs / nprocs;
    int lastsize = NumContigs - (avesize * (nprocs - 1));

    std::vector<int> recvcounts(nprocs, avesize);
    recvcounts.back() = lastsize;

    int mysize = (myrank != nprocs - 1)? avesize : lastsize;

    std::vector<int64_t> FillVecCC(mysize);

    MPI_Reduce_scatter(LocalCCSizes.data(), FillVecCC.data(), recvcounts.data(), MPI_INT64_T, MPI_SUM, World);

    FullyDistVec<int64_t, int64_t> CCSizes(FillVecCC, Assignments.getcommgrid());

    return CCSizes;
}

std::vector<std::tuple<int64_t, int64_t>> GetFilteredContigSizes(const FullyDistVec<int64_t, int64_t>& ContigSizes, int64_t MinimumContigSize)
{
    int nprocs, myrank;

    MPI_Comm World = ContigSizes.getcommgrid()->GetWorld();
    MPI_Comm_size(World, &nprocs);
    MPI_Comm_rank(World, &myrank);

    std::vector<std::tuple<int64_t, int64_t>> sendbuf;
    std::vector<int64_t> locvec = ContigSizes.GetLocVec();

    int64_t lengthuntil = ContigSizes.LengthUntil();

    for (int64_t i = 0; i < locvec.size(); ++i)
        if (locvec[i] >= MinimumContigSize)
            sendbuf.push_back(std::make_tuple(i+lengthuntil, locvec[i]));

    std::vector<int> recvcounts(nprocs);     
    std::vector<int> displs(nprocs, 0);

    recvcounts[myrank] = sendbuf.size();

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, recvcounts.data(), 1, MPI_INT, World);

    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    int64_t totalsize = std::accumulate(recvcounts.begin(), recvcounts.end(), static_cast<int64_t>(0));

    std::vector<std::tuple<int64_t, int64_t>> recvbuf(totalsize);

    MPI_Allgatherv(sendbuf.data(), recvcounts[myrank], MPIType<std::tuple<int64_t, int64_t>>(), recvbuf.data(), recvcounts.data(), displs.data(), MPIType<std::tuple<int64_t, int64_t>>(), World);

    std::sort(recvbuf.begin(), recvbuf.end(),
              [](std::tuple<int64_t, int64_t> a,
                 std::tuple<int64_t, int64_t> b) { return (std::get<1>(a) > std::get<1>(b)); });

    return recvbuf;
}

void GreedyNumberPartitioning(std::vector<std::tuple<int64_t, int64_t>> FilteredContigSizes, std::vector<int64_t>& Contig2ProcAssignments, std::vector<int64_t>& idmap)
{
    int nprocs, myrank;

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    std::vector<int64_t> sums(nprocs, 0);
    std::vector<std::vector<int64_t>> partitions(nprocs);

    for (int i = 0; i < FilteredContigSizes.size(); ++i) {
        int where = std::distance(sums.begin(), std::min_element(sums.begin(), sums.end()));
        sums[where] += std::get<1>(FilteredContigSizes[i]);
        idmap[i] = std::get<0>(FilteredContigSizes[i]);
        partitions[where].push_back(i);
    }

    for (int i = 0; i < nprocs; ++i) {
        std::vector<int64_t> mypartition = partitions[i];
        for (int j = 0; j < mypartition.size(); ++j) {
            Contig2ProcAssignments[mypartition[j]] = i;
        }
    }
}

FullyDistVec<int64_t, int64_t> GetReadAssignments
(
    FullyDistVec<int64_t, int64_t>& ContigAssignments,
    FullyDistVec<int64_t, int64_t>& ContigSizes
)
{
    int myrank, nprocs;
    MPI_Comm World = ContigAssignments.getcommgrid()->GetWorld();

    MPI_Comm_rank(World, &myrank);
    MPI_Comm_size(World, &nprocs);

    std::vector<std::tuple<int64_t, int64_t>> FilteredContigSizes = GetFilteredContigSizes(ContigSizes, 3);
    
    int64_t NumContigs = FilteredContigSizes.size();
    
    std::vector<int64_t> Contig2ProcAssignments(NumContigs);
    std::vector<int64_t> idmap(NumContigs);

    if (!myrank)
        GreedyNumberPartitioning(FilteredContigSizes, Contig2ProcAssignments, idmap);

    MPI_Bcast(Contig2ProcAssignments.data(), NumContigs, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(idmap.data(), NumContigs, MPI_INT64_T, 0, MPI_COMM_WORLD);

    std::vector<int64_t> LocalRead2ContigAssignments = ContigAssignments.GetLocVec();
    std::vector<int64_t> LocalRead2ProcAssignments(LocalRead2ContigAssignments.size(), -1);

    std::unordered_map<int64_t, int64_t> r2c;

    for (int i = 0; i < idmap.size(); ++i)
        r2c[idmap[i]] = i;

    int lengthuntil = ContigAssignments.LengthUntil();

    for (int i = 0; i < LocalRead2ContigAssignments.size(); ++i) {
        int64_t ContigID = LocalRead2ContigAssignments[i];
        if (r2c.find(ContigID) != r2c.end())
            LocalRead2ProcAssignments[i] = Contig2ProcAssignments[r2c[ContigID]];
    }

    FullyDistVec<int64_t, int64_t> Read2ProcAssignments(LocalRead2ProcAssignments, ContigAssignments.getcommgrid());
    return Read2ProcAssignments;
}


std::vector<std::vector<std::tuple<int64_t, int64_t, int64_t>>> LocalAssembly(SpDCCols<int64_t, ReadOverlap>& ContigChains, std::vector<int64_t>& LocalIdxs, int myrank)
{
    std::vector<std::vector<std::tuple<int64_t, int64_t, int64_t>>> ContigCoords;

    int64_t numreads = ContigChains.getnrow();

    assert(numreads == ContigChains.getncol());
    assert(numreads = LocalIdxs.size());

    std::vector<int64_t> roots;
    std::unordered_set<int64_t> used_roots;

    for (auto colit = ContigChains.begcol(); colit != ContigChains.endcol(); ++colit) {
        int64_t deg = 0;
        for (auto nzit = ContigChains.begnz(colit); nzit != ContigChains.endnz(colit); ++nzit, ++deg);
        if (deg==1)
            roots.push_back(colit.colid());
    }

    //for (auto itr = roots.begin(); itr != roots.end(); ++itr) {
    //    int64_t root = *itr;
    //    
    //}

    auto dcsc = ContigChains.GetDCSC();

    return ContigCoords;
}

//std::vector<std::string> LocalAssembly(SpDCCols<int64_t, ReadOverlap>& ContigChains)
//{
//    std::vector<std::string> ContigSet;
//}

#endif 
