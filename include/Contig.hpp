// Created by Giulia Guidi on 04/02/21.

// Edited by Gabriel Raulet on 12/20/21.

#ifndef CONTIG_HPP
#define CONTIG_HPP

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "CC.h"

using namespace combblas;

FullyDistVec<int64_t, int64_t> GetReadLengths(std::shared_ptr<DistributedFastaData> dfd, std::shared_ptr<CommGrid> commgrid)
{
    int myrank = commgrid->GetRank();
    int myrowrank = commgrid->GetRankInProcRow();
    int mycolrank = commgrid->GetRankInProcCol();
    int gridrows = commgrid->GetGridRows();
    int gridcols = commgrid->GetGridCols();
    int gridsize = gridrows * gridcols;

    uint64_t numreads = dfd->global_count();
    uint64_t rowreads = (myrowrank==gridrows-1)? (numreads-(gridrows-1)*(numreads/gridrows)) : (numreads/gridrows);
    uint64_t myreads  = (mycolrank==gridcols-1)? (rowreads-(gridcols-1)*(rowreads/gridcols)) : (rowreads/gridcols);

    uint64_t startidx = (rowreads/gridcols) * mycolrank;

    std::vector<int64_t> local_readlens(myreads, 0);

    for (int i = 0; i < myreads; ++i)
        local_readlens[i] = length(*(dfd->col_seq(i+startidx)));

    FullyDistVec<int64_t, int64_t> ReadLengths(local_readlens, commgrid);
    return ReadLengths;
}


FullyDistVec<int64_t, int64_t> GetContigAssignments (
    const SpParMat<int64_t, dibella::CommonKmers, SpDCCols<int64_t, dibella::CommonKmers>>& OverlapGraph,
    const FullyDistVec<int64_t, int64_t>& ReadLengths,
    FullyDistVec<int64_t, int64_t>& Branches,
    FullyDistVec<int64_t, int64_t>& Roots,
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

    Roots = degs2.Find(bind2nd(std::equal_to<int64_t>(), 1));

    FullyDistVec<int64_t, int64_t> vCC = CC(A, NumContigs);

    return vCC;
}

FullyDistVec<int64_t, int64_t> GetContigSizes (
    const SpParMat<int64_t, dibella::CommonKmers, SpDCCols<int64_t, dibella::CommonKmers>>& OverlapGraph,
    const FullyDistVec<int64_t, int64_t>& ReadLengths,
    const FullyDistVec<int64_t, int64_t>& Assignments,
    const int64_t& NumContigs
)
{
    FullyDistVec<int64_t, int64_t> overlapLens(OverlapGraph.getcommgrid());
    FullyDistVec<int64_t, int64_t> colReduce(OverlapGraph.getcommgrid());

    OverlapGraph.Reduce(overlapLens, Row, std::plus<int64_t>(), static_cast<int64_t>(0));
    OverlapGraph.Reduce(colReduce, Column, std::plus<int64_t>(), static_cast<int64_t>(0));

    overlapLens += colReduce;

    std::vector<int64_t> LocalCCSizes(NumContigs, 0);
    std::vector<int64_t> LocalCC = Assignments.GetLocVec();

    int64_t LocalCCSize = Assignments.LocArrSize();

    std::vector<int64_t> localRSizes = ReadLengths.GetLocVec();

    for (int64_t i = 0; i < LocalCCSize; ++i)
        LocalCCSizes[LocalCC[i]]++;

    int nprocs, myrank;
    MPI_Comm World = OverlapGraph.getcommgrid()->GetWorld();
    MPI_Comm_size(World, &nprocs);
    MPI_Comm_rank(World, &myrank);

    int avesize = NumContigs / nprocs;
    int lastsize = NumContigs - (avesize * (nprocs - 1));

    std::vector<int> recvcounts(nprocs, avesize);
    recvcounts.back() = lastsize;

    int mysize = (myrank != nprocs - 1)? avesize : lastsize;

    std::vector<int64_t> FillVecCC(mysize);

    MPI_Reduce_scatter(LocalCCSizes.data(), FillVecCC.data(), recvcounts.data(), MPI_INT64_T, MPI_SUM, World);

    FullyDistVec<int64_t, int64_t> CCSizes(FillVecCC, OverlapGraph.getcommgrid());

    return CCSizes;
}

#endif



