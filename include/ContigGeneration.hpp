#ifndef CONTIG_GENERATION_HPP_
#define CONTIG_GENERATION_HPP_

#include <cmath>
#include <map>
#include <fstream>
#include <cassert>
#include <limits>
#include <cstring>

#include "TraceUtils.hpp"
#include "ReadOverlap.hpp"
#include "Utils.hpp"
#include "CC.h"

using namespace combblas;

typedef int64_t IType; /* index type used in this file */

template <typename NT>
struct CSpMat
{
    typedef SpCCols<IType, NT> CCols;
    typedef SpDCCols<IType, NT> DCCols;
    typedef SpParMat<IType, NT, DCCols> MPI_DCCols;
};

typedef CSpMat<ReadOverlap>::MPI_DCCols DistStringGraph;
typedef CSpMat<ReadOverlap>::CCols LocStringGraph;
typedef FullyDistVec<IType, IType> DistAssignmentVec;

struct distinfo { MPI_Comm world; int myrank, nprocs; distinfo(std::shared_ptr<CommGrid> commgrid); };
distinfo::distinfo(std::shared_ptr<CommGrid> commgrid) : world(commgrid->GetWorld()), myrank(commgrid->GetRank()), nprocs(commgrid->GetSize()) {}

IType
GetRead2Contigs(DistStringGraph& G, DistAssignmentVec& Read2Contigs, distinfo& di);

DistAssignmentVec
GetContigSizes(const DistAssignmentVec& Read2Contigs, const IType NumContigs, distinfo& di);

DistAssignmentVec
GetRead2Procs(DistAssignmentVec& Read2Contigs, DistAssignmentVec& ContigSizes, distinfo& di);

std::vector<std::tuple<IType,IType>>
GetAllContigSizesSorted(DistAssignmentVec& ContigSizes, IType minsize, distinfo& di);

/* @func CreateContig   Assemble contigs from distributed string graph.
 *
 * @param G             combblas distributed string graph.
 * @param dfd           distributed fasta data object.
 *
 * @return vector<string> of contigs.
 */
std::vector<std::string>
CreateContig(DistStringGraph& G, std::shared_ptr<DistributedFastaData> dfd, std::string& myoutput, TraceUtils tu)
{
    float balance = G.LoadImbalance();
    int64_t nnz   = G.getnnz();

    std::ostringstream outs;
    outs << "CreateContig::LoadBalance: " << balance << std::endl;
    outs << "CreateContig::nonzeros: " << nnz << std::endl;
    SpParHelper::Print(outs.str());

    distinfo di(G.getcommgrid());

    IType NumContigs;
    DistAssignmentVec Read2Contigs;
    DistAssignmentVec ContigSizes;

    NumContigs  = GetRead2Contigs(G, Read2Contigs, di);
    ContigSizes = GetContigSizes(Read2Contigs, NumContigs, di);

    IType NumUsedContigs;
    std::vector<std::tuple<IType, IType>> AllContigSizesSorted;
    DistAssignmentVec Read2Procs;

    AllContigSizesSorted = GetAllContigSizesSorted(ContigSizes, 3, di);
    // Read2Procs  = GetRead2Procs(Read2Contigs, ContigSizes, di);

    std::vector<IType> LocalContigReadIdxs;

    LocStringGraph ContigChains = G.InducedSubgraphs2Procs(Read2Procs, LocalContigReadIdxs);

    ContigChains.Transpose();
}

/* @func GetRead2Contigs       determines which reads are in which contigs.
 *
 * @param G                    distributed string graph.
 * @Param Read2Contigs   [ref] computed assignments vector.
 *
 * @return number of contigs */
IType GetRead2Contigs(DistStringGraph& G, DistAssignmentVec& Read2Contigs, distinfo& di)
{
    CSpMat<IType>::MPI_DCCols A = G;
    CSpMat<bool>::MPI_DCCols D1, D2;
    DistAssignmentVec Branches;

    D1 = A;
    DistAssignmentVec degs1(D1.getcommgrid());

    D1.Reduce(degs1, Row, std::plus<IType>(), static_cast<IType>(0));

    Branches = degs1.FindInds(bind2nd(std::greater<IType>(), 2));

    A.PruneFull(Branches, Branches);

    D2 = A;
    DistAssignmentVec degs2(D2.getcommgrid());

    D2.Reduce(degs2, Row, std::plus<IType>(), static_cast<IType>(0));

    IType NumContigs;
    Read2Contigs = CC(A, NumContigs);

    return NumContigs;
}

/* @func GetContigSizes calculates the number of reads in each contig.
 *
 * @param Read2Contigs  read-to-contig assignments.
 * @param NumContigs    total number of contigs.
 *
 * @return distributed vector of contig sizes */
DistAssignmentVec GetContigSizes(const DistAssignmentVec& Read2Contigs, const IType NumContigs, distinfo& di)
{
    std::vector<IType> LocalCCSizes(NumContigs, 0);
    std::vector<IType> LocalCC = Read2Contigs.GetLocVec();

    IType LocalCCSize = Read2Contigs.LocArrSize();

    for (IType i = 0; i < LocalCCSize; ++i)
        LocalCCSizes[LocalCC[i]]++;

    int avesize = NumContigs / di.nprocs;
    int lastsize = NumContigs - (avesize * (di.nprocs - 1));

    std::vector<int> recvcounts(di.nprocs, avesize);
    recvcounts.back() = lastsize;

    int mysize = (di.myrank != di.nprocs - 1)? avesize : lastsize;

    std::vector<IType> FillVecCC(mysize);

    MPI_Reduce_scatter(LocalCCSizes.data(), FillVecCC.data(), recvcounts.data(), MPI_INT64_T, MPI_SUM, di.world);

    return DistAssignmentVec(FillVecCC, Read2Contigs.getcommgrid());
}

std::vector<std::tuple<IType,IType>> GetAllContigSizesSorted(DistAssignmentVec& ContigSizes, IType minsize, distinfo& di)
{
    std::vector<std::tuple<IType, IType>> sendbuf;

    std::vector<IType> locsizes = ContigSizes.GetLocVec();
    IType lengthuntil = ContigSizes.LengthUntil();

    for (IType i = 0; i < locsizes.size(); ++i)
        if (locsizes[i] >= minsize)
            sendbuf.push_back(std::make_tuple(i+lengthuntil, locsizes[i]));

    std::vector<int> recvcounts(di.nprocs);
    std::vector<int> displs(di.nprocs, 0);

    recvcounts[di.myrank] = sendbuf.size();

    MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, recvcounts.data(), 1, MPI_INT, di.world);

    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, displs.begin()+1);

    IType totalsize = std::accumulate(recvcounts.begin(), recvcounts.end(), static_cast<IType>(0));

    /* assert */

    std::vector<std::tuple<IType, IType>> result(totalsize);
    MPI_Allgatherv(sendbuf.data(), recvcounts[di.myrank], MPIType<std::tuple<IType,IType>>(), result.data(), recvcounts.data(), displs.data(), MPIType<std::tuple<IType,IType>>(), di.world);

    std::sort(result.begin(), result.end(),
              [](std::tuple<IType, IType> a,
                 std::tuple<IType, IType> b) { return (std::get<1>(a) > std::get<1>(b)); });

    return result;
}

DistAssignmentVec GetRead2Procs(DistAssignmentVec& Read2Contigs, DistAssignmentVec& ContigSizes, distinfo& di)
{

}


#endif
