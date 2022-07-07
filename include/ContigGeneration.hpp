#ifndef CONTIG_GENERATION_HPP_
#define CONTIG_GENERATION_HPP_

#include <cmath>
#include <map>
#include <fstream>
#include <cassert>
#include <limits>
#include <cstring>
#include <algorithm>

#include "TraceUtils.hpp"
#include "ReadOverlap.hpp"
#include "Utils.hpp"
#include "CC.h"

using namespace combblas;

typedef int64_t IType; /* index type used in this file */

constexpr MPI_Count max_int = std::numeric_limits<int>::max();

struct DistReadInfo
{
public:
    std::shared_ptr<CommGrid> commgrid;
    MPI_Comm world;
    int myrank, nprocs;

    FastaData *lfd;
    std::vector<IType> offsets, numreads;
    IType cached_left_val, cached_right_val;
    int cached_index;

    DistReadInfo(std::shared_ptr<CommGrid> commgrid, FastaData *lfd) : commgrid(commgrid), world(commgrid->GetWorld()), myrank(commgrid->GetRank()), nprocs(commgrid->GetSize()), lfd(lfd)
    {
        IType nlocreads = lfd->local_count();
        IType vecoffset = 0;
        MPI_Exscan(&nlocreads, &vecoffset, 1, MPIType<IType>(), MPI_SUM, world);

        offsets.resize(nprocs);
        numreads.resize(nprocs);

        MPI_Allgather(&vecoffset, 1, MPIType<IType>(), offsets.data(), 1, MPIType<IType>(), world);
        MPI_Allgather(&nlocreads, 1, MPIType<IType>(), numreads.data(), 1, MPIType<IType>(), world);

        cached_left_val = 0;
        cached_right_val = offsets.back();
        cached_index = -1;
    }

    int GetReadOwner(IType gidx)
    {
        if (gidx < 0 || gidx >= offsets.back() + numreads.back())
            return -1;

        if (cached_index != -1 && cached_left_val <= gidx && gidx < cached_right_val)
            return cached_index;


        for (IType i = 0; i < nprocs; ++i) {
            if (offsets[i] <= gidx && gidx < offsets[i] + numreads[i]) {
                cached_left_val = offsets[i];
                cached_right_val = offsets[i] + numreads[i];
                cached_index = i;
                break;
            }
        }

        return cached_index;
    }
};


IType
GetRead2Contigs(SpParMat<IType,ReadOverlap,SpDCCols<IType,ReadOverlap>>& G,
                FullyDistVec<IType,IType>&                               Read2Contigs,
                DistReadInfo&                                            di,
                TraceUtils&                                              tu);

FullyDistVec<IType,IType>
GetContigSizes(const FullyDistVec<IType,IType>& Read2Contigs,
               const IType                      NumContigs,
               DistReadInfo&                    di);

std::vector<std::tuple<IType,IType>>
GetAllContigSizesSorted(FullyDistVec<IType,IType>& ContigSizes,
                        IType&                     NumUsedContigs,
                        IType                      minsize,
                        DistReadInfo&              di,
                        TraceUtils&                tu);

std::vector<IType>
GetLocalRead2Procs(FullyDistVec<IType,IType>&            Read2Contigs,
                   std::vector<std::tuple<IType,IType>>& AllContigSizesSorted,
                   const IType                           NumUsedContigs,
                   DistReadInfo&                         di,
                   TraceUtils&                           tu);

const char *
ReadExchange(std::vector<IType>&                                   LocalRead2Procs,
             std::unordered_map<IType, std::tuple<IType, ushort>>& charbuf_info,
             DistReadInfo&                                         di,
             TraceUtils&                                           tu);

std::vector<std::string>
LocalAssembly(SpCCols<IType,ReadOverlap>&                          ContigChains,
              std::vector<IType>&                                  LocalContigReadIdxs,
              const char*                                          charbuf,
              std::unordered_map<IType, std::tuple<IType, ushort>> charbuf_info,
              DistReadInfo&                                        di);

std::vector<IType>
ImposeMyReadDistribution(FullyDistVec<IType,IType>& assignments, DistReadInfo& di);

void
AppendContig(std::string& contig, const char *buf, IType offset, ushort len, IType start, IType end);

int
MPI_Alltoallv_str(const char *sendbuf, const std::vector<IType>& sendcounts, const std::vector<IType>& sdispls,
                        char *recvbuf, const std::vector<IType>& recvcounts, const std::vector<IType>& rdispls, MPI_Comm comm);

/* @func CreateContig   Assemble contigs from distributed string graph.
 *
 * @param G             combblas distributed string graph.
 * @param dfd           distributed fasta data object.
 *
 * @return vector<string> of contigs.
 */
std::vector<std::string>
CreateContig(SpParMat<IType,ReadOverlap,SpDCCols<IType,ReadOverlap>>& G, std::shared_ptr<DistributedFastaData> dfd, std::string& myoutput, const std::shared_ptr<TimePod>& tp, TraceUtils& tu)
{
    float balance = G.LoadImbalance();
    int64_t nnz   = G.getnnz();

    std::ostringstream outs;
    outs << "CreateContig :: string graph has " << nnz << " nonzeros" << std::endl;
    tu.print_str(outs.str());
    outs.str("");

    DistReadInfo di(G.getcommgrid(), dfd->lfd());

    IType NumContigs;
    FullyDistVec<IType,IType> Read2Contigs;
    FullyDistVec<IType,IType> ContigSizes;

    tp->times["StartCreateContig:GetRead2Contigs()"] = std::chrono::system_clock::now();
    NumContigs  = GetRead2Contigs(G, Read2Contigs, di, tu);
    tp->times["EndCreateContig:GetRead2Contigs()"] = std::chrono::system_clock::now();
    tu.print_str("CreateContig :: after GetRead2Contigs\n");

    Read2Contigs.ParallelWrite("contig-assignment", true);

    tp->times["StartCreateContig:GetRead2ProcAssignments()"] = std::chrono::system_clock::now();
    ContigSizes = GetContigSizes(Read2Contigs, NumContigs, di);

    IType NumUsedContigs;
    std::vector<std::tuple<IType, IType>> AllContigSizesSorted;
    std::vector<IType> LocalRead2Procs;
    std::vector<IType> AllContig2Procs;

    AllContigSizesSorted = GetAllContigSizesSorted(ContigSizes, NumUsedContigs, 2, di, tu);

    if (!di.myrank)
    {
        std::cout << std::endl;
        std::cout << "AllContigSizesSorted:" << std::endl;
        for (auto itr = AllContigSizesSorted.begin(); itr != AllContigSizesSorted.end(); ++itr)
        {
            std::cout << std::get<0>(*itr) << "\t" << std::get<1>(*itr) << std::endl;
        }
        std::cout << std::endl;
    }

    LocalRead2Procs = GetLocalRead2Procs(Read2Contigs, AllContigSizesSorted, NumUsedContigs, di, tu);

    IType max_contig_size = ContigSizes.Reduce(combblas::maximum<IType>(), static_cast<IType>(0));
    outs << "CreateContig::NumContigs: " << NumUsedContigs << std::endl;
    outs << "CreateContig::MaxContigSize: " << max_contig_size << std::endl;
    tu.print_str(outs.str());
    outs.str("");

    FullyDistVec<IType,IType> Read2Procs(LocalRead2Procs, G.getcommgrid());
    tp->times["EndCreateContig:GetRead2ProcAssignments()"] = std::chrono::system_clock::now();
    tu.print_str("CreateContig :: Created distributed assignments vector\n");

    std::vector<IType> LocalContigReadIdxs;

    tp->times["StartCreateContig:InducedSubgraphs2Procs()"] = std::chrono::system_clock::now();
    SpDCCols<IType,ReadOverlap> ContigChainsDCC = G.InducedSubgraphs2Procs(Read2Procs, LocalContigReadIdxs);
    tp->times["EndCreateContig:InducedSubgraphs2Procs()"] = std::chrono::system_clock::now();
    tu.print_str("CreateContig :: after InducedSubgraphs2Procs\n");

    tp->times["StartCreateContig:BuildContigChains()"] = std::chrono::system_clock::now();
    SpCCols<IType,ReadOverlap> ContigChains(ContigChainsDCC);
    ContigChains.Transpose();
    tp->times["EndCreateContig:BuildContigChains()"] = std::chrono::system_clock::now();

    std::unordered_map<IType, std::tuple<IType, ushort>> charbuf_info;
    tp->times["StartCreateContig:ReadExchange()"] = std::chrono::system_clock::now();
    const char *charbuf = ReadExchange(LocalRead2Procs, charbuf_info, di, tu);
    tp->times["EndCreateContig:ReadExchange()"] = std::chrono::system_clock::now();
    tu.print_str("CreateContig :: after ReadExchange\n");

    tp->times["StartCreateContig:LocalAssembly()"] = std::chrono::system_clock::now();
    double duration = MPI_Wtime();
    std::vector<std::string> contigs = LocalAssembly(ContigChains, LocalContigReadIdxs, charbuf, charbuf_info, di);
    duration = MPI_Wtime() - duration;
    tp->times["EndCreateContig:LocalAssembly()"] = std::chrono::system_clock::now();

    double maxtime;
    MPI_Barrier(di.world);
    MPI_Reduce(&duration, &maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, di.world);

    if (!di.myrank)
    {
        std::cout << "LocalAssembly() time = " << maxtime << std::endl;
    }

    tu.print_str("CreateContig :: after LocalAssembly\n");
    delete [] charbuf;

    return contigs;
}

struct KTipsSR
{
    static IType id() { return static_cast<IType>(0); }
    static bool returnedSAID() { return false; }
    static MPI_Op mpi_op() { return MPI_LOR; }
    static IType add(const IType& arg1, const IType& arg2) { return (arg1 || arg2); }
    static IType multiply(const IType& arg1, const IType& arg2) { return (arg1 && arg2); }
    static void axpy(IType a, const IType& x, IType& y) { y = add(y, multiply(a, x)); }
};

template <class IT, class NT, class DER>
FullyDistVec<IT,IT> LastNzRowIdxPerCol(const SpParMat<IT,NT,DER>& A)
{
    std::shared_ptr<CommGrid> grid = A.getcommgrid();
    int myrank = grid->GetRank();
    int myproccol = grid->GetRankInProcRow();
    int myprocrow = grid->GetRankInProcCol();

    MPI_Comm ColWorld = grid->GetColWorld();

    IT total_rows = A.getnrow();
    IT total_cols = A.getncol();

    int procrows = grid->GetGridRows();
    int proccols = grid->GetGridCols();

    IT rows_perproc = total_rows / procrows;
    IT cols_perproc = total_cols / proccols;

    IT row_offset = myprocrow * rows_perproc;
    IT col_offset = myproccol * cols_perproc;

    DER *spSeq = A.seqptr();

    IT localcols = spSeq->getncol();
    std::vector<IT> local_colidx(localcols, static_cast<IT>(-1));

    for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
    {
        auto nzit = spSeq->begnz(colit);
        if (nzit != spSeq->endnz(colit))
            local_colidx[colit.colid()] = nzit.rowid() + row_offset;
    }

    MPI_Allreduce(MPI_IN_PLACE, local_colidx.data(), static_cast<int>(localcols), MPIType<IT>(), MPI_MAX, ColWorld);

    std::vector<IT> fillarr;

    if (!myprocrow)
        for (auto itr = local_colidx.begin(); itr != local_colidx.end(); ++itr)
            fillarr.push_back(*itr);

    return FullyDistVec<IT,IT>(fillarr, grid);
}

IType KTipsRemoval(SpParMat<IType,IType,SpDCCols<IType,IType>>& A, const FullyDistVec<IType,IType>& degrees, const IType l, TraceUtils& tu)
{
    FullyDistSpVec<IType,IType> R = degrees.Find(std::bind2nd(std::equal_to<IType>(), static_cast<IType>(1)));

    FullyDistVec<IType,IType> *ri = new FullyDistVec<IType,IType>(A.getcommgrid());
    FullyDistVec<IType,IType> *ci = new FullyDistVec<IType,IType>(A.getcommgrid());

    *ri = R.FindInds([](IType arg1) { return true; });
    ci->iota(R.getnnz(), static_cast<IType>(0));

    SpParMat<IType,IType,SpDCCols<IType,IType>> F0(A.getnrow(), R.getnnz(), *ri, *ci, static_cast<IType>(1), false);

    delete ri;
    delete ci;

    SpParMat<IType,IType,SpDCCols<IType,IType>> F1 = PSpGEMM<KTipsSR>(A, F0);
    SpParMat<IType,IType,SpDCCols<IType,IType>> V = F0;
    V += F1;

    FullyDistVec<IType,IType> TipSources(A.getcommgrid(), F0.getncol(), static_cast<IType>(-1));
    FullyDistVec<IType,IType> TipDests(A.getcommgrid(), F0.getncol(), static_cast<IType>(-1));

    for (IType k = 0; k < l; ++k)
    {
        SpParMat<IType,IType,SpDCCols<IType,IType>> F2 = PSpGEMM<KTipsSR>(A, F1);
        F2.SetDifference(V);
        V += F2;

        FullyDistVec<IType,IType> Ns = F2.Reduce(Column, std::plus<IType>(), static_cast<IType>(0));

        FullyDistSpVec<IType,IType> Tc = Ns.Find(std::bind2nd(std::greater_equal<IType>(), static_cast<IType>(2)));
        FullyDistSpVec<IType,IType> Td = Ns.Find(std::bind2nd(std::not_equal_to<IType>(), static_cast<IType>(1)));

        FullyDistVec<IType, IType> C0 = LastNzRowIdxPerCol(F0);
        FullyDistVec<IType, IType> C1 = LastNzRowIdxPerCol(F1);

        FullyDistSpVec<IType,IType> kSources = C0.GGet(Tc, [](const IType arg1, const IType arg2) { return arg2; }, static_cast<IType>(-1));
        FullyDistSpVec<IType,IType> kDests = C1.GGet(Tc, [](const IType arg1, const IType arg2) { return arg2; }, static_cast<IType>(-1));

        TipSources.Set(kSources);
        TipDests.Set(kDests);

        F1.PruneColumnByIndex(Td);
        F2.PruneColumnByIndex(Td);

        F0 = F1;
        F1 = F2;
    }

    FullyDistVec<IType,IType> where = TipSources.FindInds([](int a) { return a != -1; });
    IType num_ktip_edges = where.TotalLength();
    std::ostringstream oss;
    oss << "KTipsRemoval :: Found " << num_ktip_edges << " k-tip edges\n";
    tu.print_str(oss.str());

    A.Prune(TipSources, TipDests);
    A.Prune(TipDests, TipSources);

    return num_ktip_edges;
}

/* @func GetRead2Contigs       determines which reads belong to which which contigs.
 *
 * @param G                    distributed string graph.
 * @Param Read2Contigs   [ref] computed assignments vector.
 *
 * @description
 * Contig branching points are found by CombBLAS row reduction (to find vertex degrees),
 * and then selecting rows with more than 2 nonzeros. These are temporarily zeroed from
 * the matrix so that connected components can be run on the graph with maximum degree of 2,
 * which gives us a distributed vector assigning each vertex to a contig.
 *
 * @return number of contigs */
IType GetRead2Contigs(SpParMat<IType,ReadOverlap,SpDCCols<IType,ReadOverlap>>& G, FullyDistVec<IType,IType>& Read2Contigs, DistReadInfo& di, TraceUtils& tu)
{
    std::stringstream iss;
    SpParMat<IType,IType,SpDCCols<IType,IType>> A = G;
    tu.print_str("GetRead2Contigs :: Created IType copy of string graph for vertex degree computations\n");

    SpParMat<IType,bool,SpDCCols<IType,bool>> D1, D2;
    FullyDistVec<IType,IType> Branches;

    FullyDistVec<IType,IType> degs1(A.getcommgrid());
    FullyDistVec<IType,IType> degs2(A.getcommgrid());

    tu.print_str("GetRead2Contigs :: Calculated vertex degrees\n");

    //A.ParallelWriteMM("overlap-graph.mtx", true);

    ////IType ktip_edges_removed;
    ////do
    ////{
    ////    D1 = A;
    ////    D1.Reduce(degs1, Row, std::plus<IType>(), static_cast<IType>(0));
    ////    ktip_edges_removed = KTipsRemoval(A, degs1, 15, tu);
    ////} while (ktip_edges_removed > 0);

    //////A.ParallelWriteMM("overlap-graph-trimmed.mtx", true);

    ////tu.print_str("GetRead2Contigs :: Removed k-tips\n");

    D2 = A;
    D2.Reduce(degs2, Row, std::plus<IType>(), static_cast<IType>(0));
    tu.print_str("GetRead2Contigs :: Recalculate vertex degrees after k-tips removed\n");

    Branches = degs2.FindInds(bind2nd(std::greater<IType>(), 2));
    IType numbranches = Branches.TotalLength();
    iss << "GetRead2Contigs :: Found " << numbranches << " branching points\n";
    tu.print_str(iss.str());

    A.PruneFull(Branches, Branches);
    tu.print_str("GetRead2Contigs :: Pruned branching points\n");

    IType NumContigs;
    Read2Contigs = CC(A, NumContigs);
    iss.str("");
    iss << "GetRead2Contigs :: Found " << NumContigs << " connected components on pruned graph\n";
    tu.print_str(iss.str());

    return NumContigs;
}

/* @func GetContigSizes     calculates the number of reads within each contig.
 *
 * @param Read2Contigs      read-to-contig assignments.
 * @param NumContigs        total number of contigs.
 *
 * @return distributed vector of contig sizes */
FullyDistVec<IType,IType> GetContigSizes(const FullyDistVec<IType,IType>& Read2Contigs, const IType NumContigs, DistReadInfo& di)
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

    MPI_Reduce_scatter(LocalCCSizes.data(), FillVecCC.data(), recvcounts.data(), MPIType<IType>(), MPI_SUM, di.world);

    return FullyDistVec<IType,IType>(FillVecCC, Read2Contigs.getcommgrid());
}
/* @func GetAllContigSizesSorted
 *
 * @param ContigSizes
 * @param NumUsedContigs
 * @param minsize
 *
 * @description
 * Using the distributed vector of contig sizes as input, find all contigs with
 * >= @minsize reads, and All-gather them so that each processor has a vector
 * of (contig-id (global), contig-size) tuples. Then each processor sorts their
 * vector by contig-size and returns the result.
 *
 * @return */
std::vector<std::tuple<IType,IType>> GetAllContigSizesSorted(FullyDistVec<IType,IType>& ContigSizes, IType& NumUsedContigs, IType minsize, DistReadInfo& di, TraceUtils& tu)
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

    std::ostringstream outs;
    outs << "GetAllContigSizesSorted :: Found " << NumUsedContigs << " connected components with size >= 2\n";
    tu.print_str(outs.str());

    std::vector<std::tuple<IType, IType>> result(NumUsedContigs);
    MPI_Allgatherv(sendbuf.data(), recvcounts[di.myrank], MPIType<std::tuple<IType,IType>>(), result.data(), recvcounts.data(), displs.data(), MPIType<std::tuple<IType,IType>>(), di.world);
    tu.print_str("GetAllContigSizesSorted :: All gathered contig ids and sizes\n");

    std::sort(result.begin(), result.end(),
              [](std::tuple<IType, IType> a,
                 std::tuple<IType, IType> b) { return (std::get<1>(a) > std::get<1>(b)); });
    tu.print_str("GetAllContigSizesSorted :: Sorted contig ids by size\n");

    return result;
}

/* @func ImposeMyReadDistribution
 *
 * @param assignments combblas distributed assignments vector
 *
 * @description
 * Takes a combblas distributed assignments vector and returns the local section
 * of the vector redistributed according to the read distribution of diBELLA.
 *
 * @return redistributed local section of assignments vector.
 */
std::vector<IType> ImposeMyReadDistribution(FullyDistVec<IType,IType>& assignments, DistReadInfo& di)
{
    IType orig_offset = assignments.LengthUntil();
    IType orig_size = assignments.LocArrSize();

    std::vector<IType> orig_vector = assignments.GetLocVec();

    std::vector<int> sendcounts(di.nprocs, 0);
    std::vector<int> recvcounts(di.nprocs);

    for (IType i = 0; i < orig_size; ++i) {
        int owner = di.GetReadOwner(i+orig_offset);
        ++sendcounts[owner];
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, di.world);

    std::vector<int> sdispls(di.nprocs, 0);
    std::vector<int> rdispls(di.nprocs, 0);

    std::partial_sum(sendcounts.begin(), sendcounts.end()-1, sdispls.begin()+1);
    std::partial_sum(recvcounts.begin(), recvcounts.end()-1, rdispls.begin()+1);

    IType totrecv = std::accumulate(recvcounts.begin(), recvcounts.end(), static_cast<IType>(0));

    assert((totrecv < std::numeric_limits<int>::max()));

    std::vector<IType> new_vector(totrecv);

    MPI_Alltoallv(orig_vector.data(), sendcounts.data(), sdispls.data(), MPIType<IType>(), new_vector.data(), recvcounts.data(), rdispls.data(), MPIType<IType>(), di.world);

    return new_vector;
}

/* @func GetLocalRead2Procs      determine the local read to processor assignments for load
 *                               balancing contig generation.
 *
 * @param Read2Contigs           distributed read to contig id assignments.
 * @param AllContigSizesSorted   global vector of (contig id, size) tuples, sorted by size.
 * @param
 *
 * @description
 * Using the @Read2Contigs assignments and globally shared @AllContigSizesSorted vector of
 * used contig ids (and sizes), we compute the partitioning of used contigs to processor ranks using
 * greedy multiway number partitioning optimization algorithm. We call 'small'-indices the
 * indices of used contigs within the smaller @AllContigSizesSorted vector, and we call
 * 'large'-indices the corresponding global contig ids originally calculated by connected
 * components. The small indices are used for the optimization algorithm since they enable
 * us a headache-free way to run this algorithm and to broadcast the partitioning to
 * every processor (note: the optimization algorithm is only run on the root rank).
 *
 * Once the optimization algorithm is finished on the root rank, the used contig-to-processor
 * assignments are broadcast to everyone, and the mapping of the small indices to the large
 * indices is broadcast as well.
 *
 * We then need to create a hashmap which inverts the small-to-large map and gives us a
 * large-to-small map, so that we know which global contig-id we are referring to when
 * we are looking at a small index.
 *
 * At this point, we want to finally create a local vector for each processor which tells
 * us where the reads that our processor owns are supposed to end up. The idea is that since
 * we know which reads belong to which contigs, and we now know which contigs are to be sent
 * to which processors, we automatically know which reads are to be sent to which processors.
 * So we simply go through our local read-to-contig assignments and, using the globally computed
 * contig partititioning, determine where our local reads are supposed to go.
 *
 * @return local read to processor assignments */
std::vector<IType> GetLocalRead2Procs(FullyDistVec<IType,IType>& Read2Contigs, std::vector<std::tuple<IType,IType>>& AllContigSizesSorted, const IType NumUsedContigs, DistReadInfo& di, TraceUtils& tu)
{
    assert((NumUsedContigs == static_cast<IType>(AllContigSizesSorted.size())));

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

        std::cout << std::endl;
        std::cout << "SmallToLargeMap:" << std::endl;
        for (IType i = 0; i < SmallToLargeMap.size(); ++i)
        {
            std::cout << i << "\t" << SmallToLargeMap[i] << std::endl;
        }
        std::cout << std::endl;
    }

    tu.print_str("GetLocalRead2Procs :: Root process finished multiway partitioning for contig load-balancing\n");

    MPI_Bcast(AllContig2Procs.data(), NumUsedContigs, MPIType<IType>(), 0, di.world);
    MPI_Bcast(SmallToLargeMap.data(), NumUsedContigs, MPIType<IType>(), 0, di.world);

    tu.print_str("GetLocalRead2Procs :: Broadcast contig-to-processor assignments\n");

    std::unordered_map<IType,IType> LargeToSmallMap; /* maps all global 'large' contig ids to 'small' ids */
    for (auto itr = SmallToLargeMap.begin(); itr != SmallToLargeMap.end(); ++itr)
        LargeToSmallMap[*itr] = (itr - SmallToLargeMap.begin());

    IType lengthuntil = di.offsets[di.myrank];
    IType nlocreads = di.numreads[di.myrank];

    std::vector<IType> LocalRead2Contigs = ImposeMyReadDistribution(Read2Contigs, di);
    tu.print_str("GetLocalRead2Procs :: Got reconfigured read-to-processor distribution (local segment)\n");

    assert((nlocreads == static_cast<IType>(LocalRead2Contigs.size())));

    std::vector<IType> LocalRead2Procs(nlocreads, -1);

    for (auto itr = LocalRead2Contigs.begin(); itr != LocalRead2Contigs.end(); ++itr)
       if (LargeToSmallMap.find(*itr) != LargeToSmallMap.end())
           LocalRead2Procs[itr - LocalRead2Contigs.begin()] = AllContig2Procs[LargeToSmallMap[*itr]];

    tu.print_str("GetLocalRead2Procs :: Got local read-to-processor assignments\n");

    return LocalRead2Procs;
}

/* @func ReadExchange          communicate reads that each processor doesn't have access to,
                               but which needs in order to generate contigs.
 *
 * @param LocalRead2Procs      local read-to-processor assignments.
 * @param charbuf_info [ref]   hashmap of global read-indices to (charbuffer-offset, read-length) tuples.
 *
 * @description
 * Plan is to send a packed char buffer of reads to each processor and to send a
 * map @charbuf_info which has the associated offset within the char buffer and the
 * length of the read, so that read sequences can be quickly looked up within the
 * char buffer. We do this by communicating two objects: the char buffer,
 * and an array of info-tuples ordered according to how the char buffer itself is ordered.
 *
 * First step is to loop through @LocalRead2Procs vector for reads which need
 * to be sent to a different processor. For each such read, we collect the information
 * regarding where it is to be sent, where its offset is within the @lfd_buffer, and
 * how long it is. We also calculate the char buffer sendcounts in the same loop.
 * We then pack up the tuple and counts information and send them around.
 *
 * Then we allocate the char send and receive buffers because we now know how long
 * they should each be, and then we fill the char send buffer. We then send
 * the char buffers around.
 *
 * Finally, we use the tuple buffers to construct the local @charbuf_info map on
 * each process, and then we are done.
 *
 * @return char buffer pointer
 */

const char * ReadExchange(std::vector<IType>& LocalRead2Procs, std::unordered_map<IType, std::tuple<IType, ushort>>& charbuf_info, DistReadInfo& di, TraceUtils& tu)
{
    IType char_totsend = 0, char_totrecv;
    IType read_totsend = 0, read_totrecv;

    IType nlocreads = di.numreads[di.myrank];
    IType vecoffset  = di.offsets[di.myrank];

    char *char_send, *char_recv;
    std::tuple<IType, IType> *read_sendbuf, *read_recvbuf;

    std::vector<int>   read_sendcounts (di.nprocs, 0);
    std::vector<int>   read_recvcounts (di.nprocs   );
    std::vector<int>   read_sdispls    (di.nprocs, 0);
    std::vector<int>   read_rdispls    (di.nprocs, 0);
    std::vector<IType> char_sendcounts (di.nprocs, 0);
    std::vector<IType> char_recvcounts (di.nprocs   );
    std::vector<IType> char_sdispls    (di.nprocs, 0);
    std::vector<IType> char_rdispls    (di.nprocs, 0);

    assert((nlocreads == static_cast<IType>(LocalRead2Procs.size())));

    std::vector<std::vector<std::tuple<IType, IType>>> i_read_sendbuf(di.nprocs);
    std::unordered_map<IType, std::tuple<ushort, IType>> lfd_info_cache;

    /* read_sendbuf is 'vector' of tuples (global read index, read length) */

    for (IType i = 0; i < nlocreads; ++i) {
        int dest = LocalRead2Procs[i]; /* dest == -1 means that the read is not in a contig */
        if (dest >= 0 && di.myrank != dest) {
            ushort readlen;
            uint64_t start_offset, end_offset;
            di.lfd->get_sequence(i, readlen, start_offset, end_offset);
            lfd_info_cache[i] = std::make_tuple(readlen, start_offset);
            i_read_sendbuf[dest].push_back(std::make_tuple(i+vecoffset, readlen));
            char_sendcounts[dest] += readlen;
            char_totsend += readlen;
            ++read_sendcounts[dest];
            ++read_totsend;
        }
    }

    read_sendbuf = new std::tuple<IType, IType>[read_totsend];

    std::partial_sum(read_sendcounts.begin(), read_sendcounts.end()-1, read_sdispls.begin()+1);
    std::partial_sum(char_sendcounts.begin(), char_sendcounts.end()-1, char_sdispls.begin()+1);

    for (int i = 0; i < di.nprocs; ++i)
        std::copy(i_read_sendbuf[i].begin(), i_read_sendbuf[i].end(), read_sendbuf + read_sdispls[i]);

    tu.print_str("ReadExchange :: Filled read id/size send buffers\n");

    MPI_Alltoall(read_sendcounts.data(), 1, MPI_INT,     read_recvcounts.data(), 1, MPI_INT,     di.world);
    MPI_Alltoall(char_sendcounts.data(), 1, MPI_INT64_T, char_recvcounts.data(), 1, MPI_INT64_T, di.world);

    read_totrecv = std::accumulate(read_recvcounts.begin(), read_recvcounts.end(), static_cast<IType>(0));
    char_totrecv = std::accumulate(char_recvcounts.begin(), char_recvcounts.end(), static_cast<IType>(0));

    std::partial_sum(read_recvcounts.begin(), read_recvcounts.end()-1, read_rdispls.begin()+1);
    std::partial_sum(char_recvcounts.begin(), char_recvcounts.end()-1, char_rdispls.begin()+1);

    read_recvbuf = new std::tuple<IType,IType>[read_totrecv];

    assert((read_totrecv < std::numeric_limits<int>::max()));
    assert((read_totsend < std::numeric_limits<int>::max()));

    MPI_Alltoallv(read_sendbuf, read_sendcounts.data(), read_sdispls.data(), MPIType<std::tuple<IType,IType>>(),
                  read_recvbuf, read_recvcounts.data(), read_rdispls.data(), MPIType<std::tuple<IType,IType>>(), di.world);

    tu.print_str("ReadExchange :: All-to-all communicated read id/size information\n");

    char_send = new char[char_totsend+1]();
    char_recv = new char[char_totrecv+1](); /* allocated charbuf */

    ushort len;
    IType offset;
    const char *lfd_buffer = di.lfd->buffer();
    char *char_send_ptr = char_send;

    for (IType i = 0; i < read_totsend; ++i) {
        auto cached = lfd_info_cache[std::get<0>(read_sendbuf[i]) - vecoffset];
        len = std::get<0>(cached), offset = std::get<1>(cached);
        std::strncpy(char_send_ptr, lfd_buffer + offset, len);
        char_send_ptr += len;
    }

    tu.print_str("ReadExchange :: Filled char buffers\n");

    MPI_Alltoallv_str(char_send, char_sendcounts, char_sdispls, char_recv, char_recvcounts, char_rdispls, di.world);
    tu.print_str("ReadExchange :: Completed custom all-to-all using derived datatypes and Isend/Irecv calls\n");

    /* TODO: inefficient */
    std::vector<IType> char_read_offsets(read_totrecv, 0);
    for (IType i = 0; i < read_totrecv-1; ++i)
        char_read_offsets[i+1] = char_read_offsets[i] + std::get<1>(read_recvbuf[i]);

    for (IType i = 0; i < read_totrecv; ++i)
        charbuf_info[std::get<0>(read_recvbuf[i])] = std::make_tuple(char_read_offsets[i], std::get<1>(read_recvbuf[i]));

    tu.print_str("ReadExchange :: Filled out char buffer information maps for read sequence lookup\n");
    delete [] char_send;
    delete [] read_sendbuf;
    delete [] read_recvbuf;

    return char_recv;
}

/* @func LocalAssembly           Assemble contigs from locally induced subgraph and read sequence information.
 *
 * @param ContigChains           locally induced subgraph of contig chains.
 * @param LocalContigReadIdxs    mapping of local graph indices to original global read indices.
 * @param charbuf                char buffer of received reads (not including reads already here).
 * @param charbuf_info           information for char buffer lookup/reading.
 *
 * @description
 *
 * @returns vector of completed contigs */
std::vector<std::string> LocalAssembly(SpCCols<IType,ReadOverlap>& ContigChains, std::vector<IType>& LocalContigReadIdxs, const char *charbuf, std::unordered_map<IType, std::tuple<IType, ushort>> charbuf_info, DistReadInfo& di)
{
    /* local fasta buffer for sequences already on my processor */
    const char *lfd_buffer = di.lfd->buffer();

    std::vector<std::string> contigs;

    IType vecoffset = di.offsets[di.myrank];
    IType numreads = ContigChains.getnrow();
    assert((numreads = static_cast<IType>(ContigChains.getncol())));
    assert((numreads = static_cast<IType>(LocalContigReadIdxs.size())));

    auto csc = ContigChains.GetCSC();

    if (!csc)
        return std::vector<std::string>();

    assert((numreads>=2));

    bool *visited = new bool[numreads](); /* vector of visited reads */
    std::unordered_set<IType> used_roots; /* because there are two roots per contig, and we only
                                             want to use one each, record the which ones have been
                                             used so we avoid two traversals per contig */

    /* we loop through each vertex searching for roots */
    for (IType v = 0; v < csc->n; ++v) {

        /* @jc is the column pointer vector. Because roots are defined as
         * as vertices of degree 1, we check whether the column size is 1,
         * and if so, whether this root has already been traversed by a
         * previous contig */
        if (csc->jc[v+1] - csc->jc[v] != 1 || used_roots.find(v) != used_roots.end())
            continue;

        /* for each read that participates in a contig, we store the start and end indices
         * of the correct substring as well as the gobal read index */
        std::vector<std::tuple<IType, IType, IType>> contig_vector;

        IType cur = v;
        IType end, next;
        IType i1last = 0;

        ReadOverlap e;

        bool first = true;

        while (true) {
            visited[cur] = true;
            next = csc->jc[cur];
            end = csc->jc[cur+1];

            while (next < end && visited[csc->ir[next]])
                ++next;

            if (next >= end)
                break;

            e = csc->num[next];

            if (first) {
                i1last = (e.dir == 0 || e.dir == 1)? 0 : e.l[0];
                first = false;
            }

            contig_vector.push_back(std::make_tuple(i1last, e.coords[0], LocalContigReadIdxs[cur]));

            i1last = e.coords[1];
            cur = csc->ir[next];
        }

        contig_vector.push_back(std::make_tuple(i1last, (e.dir == 1 || e.dir == 3)? e.l[1] : 0, LocalContigReadIdxs[cur]));
        used_roots.insert(cur);

        std::string contig = "";

        for (IType i = 0; i < contig_vector.size(); ++i) {

            IType seq_start = std::get<0>(contig_vector[i]);
            IType seq_end   = std::get<1>(contig_vector[i]);
            IType idx       = std::get<2>(contig_vector[i]);

            auto segment_info_itr = charbuf_info.find(idx);
            uint64_t read_offset, end_offset;
            ushort readlen;

            if (segment_info_itr == charbuf_info.end()) {
                di.lfd->get_sequence(idx-vecoffset, readlen, read_offset, end_offset);
                AppendContig(contig, lfd_buffer, read_offset, readlen, seq_start, seq_end);
            } else {
                std::tuple<IType, ushort> val = segment_info_itr->second;
                read_offset = std::get<0>(val), readlen = std::get<1>(val);
                AppendContig(contig, charbuf, read_offset, readlen, seq_start, seq_end);
            }
        }
        contigs.push_back(contig);
    }

    delete [] visited;

    return contigs;
}

// 0000000000111111111122222222223333333333
// 0123456789012345678901234567890123456789
// ACCGCCGTATCGCGCGATATATATGCGCGTAGAGCGCCCA
// 9876543210987654321098765432109876543210
// 3333333333222222222211111111110000000000
//
//           00000000
//           23456789
// w[2,10] = CGCCGTAT
//
//
// w[10,2] = ATACGGCG = comp(w)[len(w)-10, len(w)-2] = comp(w)[30,38]
//           76543210
//           33333333

char comp(char c)
{
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
    }
    return '\0';
}

void AppendContig(std::string& contig, const char *buf, IType offset, ushort len, IType start, IType end)
{
    if (start < end)
        contig.append(buf + offset + start, end - start);
    else {
        std::string segment(buf + offset + end, start - end);
        std::transform(segment.cbegin(), segment.cend(), segment.begin(), comp);
        std::reverse(segment.begin(), segment.end());
        contig.append(segment);
    }
}

int MPI_Alltoallv_str(const char *sendbuf, const std::vector<IType>& sendcounts, const std::vector<IType>& sdispls,
                            char *recvbuf, const std::vector<IType>& recvcounts, const std::vector<IType>& rdispls, MPI_Comm comm)
{
    assert((sizeof(void*)==8));

    int nprocs, myrank;
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);

    MPI_Request *reqs = new MPI_Request[2*nprocs];

    for (int i = 0; i < nprocs; ++i) {
        void *buf = recvbuf + rdispls[i];
        MPI_Count recvcount = recvcounts[i];
        if (recvcount <= max_int) {
            MPI_Irecv(buf, static_cast<int>(recvcount), MPI_CHAR, i, 0, comm, &reqs[i]);
        } else {
            MPI_Datatype newtype, chunks, remains;
            MPI_Type_vector(recvcount/max_int, max_int, max_int, MPI_CHAR, &chunks);
            MPI_Type_contiguous(recvcount%max_int, MPI_CHAR, &remains);
            MPI_Aint rdisp = static_cast<MPI_Aint>((recvcount/max_int)*max_int);
            int blklens[2] = {1,1};
            MPI_Aint displs[2] = {0,rdisp};
            MPI_Datatype types[2] = {chunks,remains};
            MPI_Type_create_struct(2, blklens, displs, types, &newtype);
            MPI_Type_free(&chunks);
            MPI_Type_free(&remains);
            MPI_Type_commit(&newtype);
            MPI_Irecv(buf, 1, newtype, i, 0, comm, &reqs[i]);
            MPI_Type_free(&newtype);
        }
    }

    for (int j = myrank; j < (nprocs+myrank); ++j) {
        int i = j%nprocs;
        const void *buf = sendbuf + sdispls[i];
        MPI_Count sendcount = sendcounts[i];
        if (sendcount <= max_int) {
            MPI_Isend(buf, static_cast<int>(sendcount), MPI_CHAR, i, 0, comm, &reqs[i+nprocs]);
        } else {
            MPI_Datatype newtype, chunks, remains;
            MPI_Type_vector(sendcount/max_int, max_int, max_int, MPI_CHAR, &chunks);
            MPI_Type_contiguous(sendcount%max_int, MPI_CHAR, &remains);
            MPI_Aint rdisp = static_cast<MPI_Aint>((sendcount/max_int)*max_int);
            int blklens[2] = {1,1};
            MPI_Aint displs[2] = {0,rdisp};
            MPI_Datatype types[2] = {chunks,remains};
            MPI_Type_create_struct(2, blklens, displs, types, &newtype);
            MPI_Type_free(&chunks);
            MPI_Type_free(&remains);
            MPI_Type_commit(&newtype);
            MPI_Isend(buf, 1, newtype, i, 0, comm, &reqs[i+nprocs]);
            MPI_Type_free(&newtype);
        }
    }

    MPI_Waitall(2*nprocs, reqs, MPI_STATUSES_IGNORE);
    delete [] reqs;

    return MPI_SUCCESS;
}

#endif
