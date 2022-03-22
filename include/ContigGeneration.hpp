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

struct ReadLocations;

struct distinfo { MPI_Comm world; int myrank, nprocs; distinfo(std::shared_ptr<CommGrid> commgrid); };
distinfo::distinfo(std::shared_ptr<CommGrid> commgrid) : world(commgrid->GetWorld()), myrank(commgrid->GetRank()), nprocs(commgrid->GetSize()) {}

struct ReadLocations
{
public:
    std::vector<IType> offsets;
    std::vector<IType> numreads;

    IType cached_left_val, cached_right_val;
    int cached_index;

    ReadLocations(const IType num_locreads, const distinfo& di)
    {
        IType myoffset = 0;
        MPI_Exscan(&num_locreads, &myoffset, 1, MPIType<IType>(), MPI_SUM, di.world);

        offsets.resize(di.nprocs);
        numreads.resize(di.nprocs);

        MPI_Allgather(&myoffset, 1, MPIType<IType>(), offsets.data(), 1, MPIType<IType>(), di.world);
        MPI_Allgather(&num_locreads, 1, MPIType<IType>(), numreads.data(), 1, MPIType<IType>(), di.world);

        cached_left_val = 0;
        cached_right_val = offsets.back();
        cached_index = -1;
    }

    int GetReadOwner(IType gidx)
    {
        if (cached_left_val <= gidx && gidx < cached_right_val)
            return cached_index;

        if (gidx < 0 || gidx >= offsets.back())
            return -1;

        int left, right, mid;
        left = 0, right = offsets.size() - 1;

        while (left < right) {
            mid = (left + right) / 2;
            if (offsets[mid] <= gidx && gidx < offsets[mid] + numreads[mid])
                break;
            else if (offsets[mid] < gidx)
                left = mid + 1;
            else if (offsets[mid] > gidx)
                right = mid - 1;
            else
                break;
        }
        cached_left_val = offsets[mid];
        cached_right_val = cached_left_val + numreads[mid];
        cached_index = mid;

        return cached_index;
    }
};

IType
GetRead2Contigs(DistStringGraph& G, DistAssignmentVec& Read2Contigs, distinfo& di);

DistAssignmentVec
GetContigSizes(const DistAssignmentVec& Read2Contigs, const IType NumContigs, distinfo& di);

std::vector<std::tuple<IType,IType>>
GetAllContigSizesSorted(DistAssignmentVec& ContigSizes, IType& NumUsedContigs, IType minsize, distinfo& di);

std::vector<IType>
GetLocalRead2Procs(DistAssignmentVec& Read2Contigs, std::vector<std::tuple<IType,IType>>& AllContigSizesSorted, const IType NumUsedContigs, ReadLocations& readlocs, distinfo& di);


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
    ReadLocations readlocs(dfd->lfd()->local_count(), di);

    IType NumContigs;
    DistAssignmentVec Read2Contigs;
    DistAssignmentVec ContigSizes;

    NumContigs  = GetRead2Contigs(G, Read2Contigs, di);
    ContigSizes = GetContigSizes(Read2Contigs, NumContigs, di);

    IType NumUsedContigs;
    std::vector<std::tuple<IType, IType>> AllContigSizesSorted;
    std::vector<IType> LocalRead2Procs;
    std::vector<IType> AllContig2Procs;

    AllContigSizesSorted = GetAllContigSizesSorted(ContigSizes, NumUsedContigs, 3, di);
    LocalRead2Procs      = GetLocalRead2Procs(Read2Contigs, AllContigSizesSorted, NumUsedContigs, readlocs, di);

    std::vector<IType> LocalContigReadIdxs;
    LocStringGraph ContigChains = G.InducedSubgraphs2Procs(DistAssignmentVec(LocalRead2Procs, G.getcommgrid()), LocalContigReadIdxs);

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

std::vector<std::tuple<IType,IType>> GetAllContigSizesSorted(DistAssignmentVec& ContigSizes, IType& NumUsedContigs, IType minsize, distinfo& di)
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

    NumUsedContigs = std::accumulate(recvcounts.begin(), recvcounts.end(), static_cast<IType>(0));

    /* assert */

    std::vector<std::tuple<IType, IType>> result(NumUsedContigs);
    MPI_Allgatherv(sendbuf.data(), recvcounts[di.myrank], MPIType<std::tuple<IType,IType>>(), result.data(), recvcounts.data(), displs.data(), MPIType<std::tuple<IType,IType>>(), di.world);

    std::sort(result.begin(), result.end(),
              [](std::tuple<IType, IType> a,
                 std::tuple<IType, IType> b) { return (std::get<1>(a) > std::get<1>(b)); });

    return result;
}


/* @func GetLocalRead2Procs      determine the local read to processor assignments for load
 *                               balancing contig generation.
 *
 * @param Read2Contigs           distributed read to contig id assignments.
 * @param AllContigSizesSorted   global vector of (contig id, size) tuples, sorted by size.
 * @param NumUsedContigs         number of contigs that will be generated.
 * @param readlocs               read locations data structure for read owner query.
 *
 * @return local read to processor assignments */
std::vector<IType> GetLocalRead2Procs(DistAssignmentVec& Read2Contigs, std::vector<std::tuple<IType,IType>>& AllContigSizesSorted, const IType NumUsedContigs, ReadLocations& readlocs, distinfo& di)
{
    std::vector<IType> AllContig2Procs(NumUsedContigs); /* maps all global 'small' contig ids to their destination processor (every rank has this same copy after MPI_Bcast below */
    std::vector<IType> SmallToLargeMap(NumUsedContigs); /* maps compressed used contig ids (since NumUsedContigs << NumContigs) to distributed ContigIDs */

    /* multiway number partitioning, greedy approach */
    if (!di.myrank) {
        std::vector<IType> sums(di.nprocs, 0);
        std::vector<std::vector<IType>> partitions(di.nprocs);

        for (IType i = 0; i < AllContigSizesSorted.size(); ++i) {
            int where = std::distance(sums.begin(), std::min_element(sums.begin(), sums.end()));
            sums[where] += std::get<1>(AllContigSizesSorted[i]);
            SmallToLargeMap[i] = std::get<0>(AllContigSizesSorted[i]);
            partitions[where].push_back(i);
        }

        for (int i = 0; i < di.nprocs; ++i) {
            std::vector<IType> mypartition = partitions[i];
            for (auto itr = mypartition.begin(); itr != mypartition.end(); ++itr)
                AllContig2Procs[*itr] = i;
        }
    }

    MPI_Bcast(AllContig2Procs.data(), NumUsedContigs, MPIType<IType>(), 0, di.world);
    MPI_Bcast(SmallToLargeMap.data(), NumUsedContigs, MPIType<IType>(), 0, di.world);

    std::unordered_map<IType,IType> LargeToSmallMap; /* maps all global 'large' contig ids to 'small' ids */
    for (auto itr = SmallToLargeMap.begin(); itr != SmallToLargeMap.end(); ++itr)
        LargeToSmallMap[*itr] = (itr - SmallToLargeMap.begin());

    IType myoffset = Read2Contigs.LengthUntil();
    IType mysize = Read2Contigs.LocArrSize();

    std::vector<int> sendcounts(di.nprocs, 0);
    std::vector<int> recvcounts(di.nprocs);

    for (IType i = myoffset; i < myoffset + mysize; ++i) {
        int owner = readlocs.GetReadOwner(i);
        ++sendcounts[owner];
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, di.world);

    std::vector<int> sdispls(di.nprocs, 0);
    std::vector<int> rdispls(di.nprocs, 0);

    std::partial_sum(sendcounts.begin(), sendcounts.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, rdispls.begin()+1);

    IType totrecv = std::accumulate(recvcounts.begin(), recvcounts.end(), static_cast<IType>(0));

    std::vector<IType> LocalRead2Contigs_combblas = Read2Contigs.GetLocVec();
    std::vector<IType> LocalRead2Contigs(totrecv);

    MPI_Alltoallv(LocalRead2Contigs_combblas.data(), sendcounts.data(), sdispls.data(), MPIType<IType>(), LocalRead2Contigs.data(), recvcounts.data(), rdispls.data(), MPIType<IType>(), di.world);

    std::vector<IType> LocalRead2Procs(LocalRead2Contigs.size(), -1);

    IType lengthuntil = readlocs.offsets[di.myrank];
    for (auto itr = LocalRead2Contigs.begin(); itr != LocalRead2Contigs.end(); ++itr)
       if (LargeToSmallMap.find(*itr) != LargeToSmallMap.end())
           LocalRead2Procs[itr - LocalRead2Contigs.begin()] = AllContig2Procs[LargeToSmallMap[*itr]];

    return LocalRead2Procs;
}

#endif
