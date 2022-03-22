
#ifndef CONTIG_HPP
#define CONTIG_HPP

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

template <typename T>
void write_object_to_file(T object, std::string prefix)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    std::stringstream filename;
    filename << prefix << "_" << myrank << ".txt";
    std::ofstream filestream(filename.str());
    filestream << object << std::endl;
    filestream.close();
}

template <typename T>
void write_vector_to_file(const std::vector<T>& vec, std::string prefix)
{
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    std::stringstream filename;
    filename << prefix << "_" << myrank << ".txt";
    std::ofstream filestream(filename.str());
    for (auto itr = vec.begin(); itr != vec.end(); ++itr) {
        filestream << *itr << ", ";
    }
    filestream << std::endl;
    filestream.close();
}

FullyDistVec<int64_t, int64_t> GetContigAssignments
(
    const SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>& OverlapGraph,
    FullyDistVec<int64_t, int64_t>& Branches,
    FullyDistVec<int64_t, int64_t>& Roots,
    int64_t& NumContigs
)
{
    /* copy */
    SpParMat<int64_t, int64_t, SpDCCols<int64_t, int64_t>> A = OverlapGraph;

    /* boolean matrix used for branch finding calculation */
    SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> D1 = A;
    FullyDistVec<int64_t, int64_t> degs1(D1.getcommgrid());

    /* get vertex degrees */
    D1.Reduce(degs1, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    /* branchces are those vertices of degree 3 or greater */
    Branches = degs1.FindInds(bind2nd(std::greater<int64_t>(), 2));

    /* get rid of edges incident with all branching vertices */
    A.PruneFull(Branches, Branches);

    /* boolean matrix used for root finding calculation */
    SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> D2 = A;
    FullyDistVec<int64_t, int64_t> degs2(D2.getcommgrid());

    D2.Reduce(degs2, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    /* root vertices have degree 1 */
    Roots = degs2.FindInds(bind2nd(std::equal_to<int64_t>(), 1));

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

    assert((totalsize < std::numeric_limits<int>::max()));

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

std::vector<std::vector<std::tuple<int64_t, int64_t, int64_t>>> LocalAssembly(SpCCols<int64_t, ReadOverlap>& ContigChains, std::vector<int64_t>& LocalIdxs, int myrank)
{
    std::vector<std::vector<std::tuple<int64_t, int64_t, int64_t>>> ContigCoords;

    int64_t numreads = ContigChains.getnrow();

    assert((numreads == ContigChains.getncol()));
    assert((numreads = LocalIdxs.size()));

    bool visited[numreads] = {0};
    std::unordered_set<int64_t> used_roots;

    auto csc = ContigChains.GetCSC();

    for (int64_t v = 0; v < csc->n; ++v) {

        if (csc->jc[v+1] - csc->jc[v] != 1 || used_roots.find(v) != used_roots.end())
            continue;

        std::vector<std::tuple<int64_t, int64_t, int64_t>> contig;

        int64_t cur = v;
        int64_t end, next;
        int64_t i1last = 0;

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

            contig.push_back(std::make_tuple(i1last, e.coords[0], LocalIdxs[cur]));

            i1last = e.coords[1];
            cur = csc->ir[next];
        }

        contig.push_back(std::make_tuple(i1last, (e.dir == 1 || e.dir == 3)? e.l[1] : 0, LocalIdxs[cur]));

        used_roots.insert(cur);
        ContigCoords.push_back(contig);
    }

    return ContigCoords;
}

int64_t ReadExchange(FullyDistVec<int64_t, int64_t> Read2ProcAssignments,
                     const FastaData *lfd,
                     char **charbuf,
                     std::vector<uint64_t>& charbuf_read_idxs,
                     std::vector<uint64_t>& charbuf_read_lengths,
                     std::vector<uint64_t>& charbuf_read_offsets,
                     uint64_t& myoffset)
{
    MPI_Comm World = Read2ProcAssignments.getcommgrid()->GetWorld();

    int myrank, nprocs;
    MPI_Comm_rank(World, &myrank);
    MPI_Comm_size(World, &nprocs);

    std::vector<int64_t> read2procs = Read2ProcAssignments.GetLocVec(); /* local vector segment of distributed @Read2ProcAssignments,
                                                                           mapping local read indices to their destination processes */
    //uint64_t myoffset = Read2ProcAssignments.LengthUntil(); /* global offset of @Read2ProcAssignments */
    //uint64_t num_locreads = read2procs.size(); /* number of reads stored on local processor */

    std::vector<int64_t> read2procs_global(Read2ProcAssignments.TotalLength());

    std::vector<int> read2procs_recvcounts(nprocs);
    std::vector<int> read2procs_displs(nprocs, 0);

    int s = static_cast<int>(read2procs.size());
    MPI_Allgather(&s, 1, MPI_INT, read2procs_recvcounts.data(),  1, MPI_INT, World);

    std::partial_sum(read2procs_recvcounts.begin(), read2procs_recvcounts.end()-1, read2procs_displs.begin()+1);

    MPI_Allgatherv(read2procs.data(), s, MPI_INT64_T, read2procs_global.data(), read2procs_recvcounts.data(), read2procs_displs.data(), MPI_INT64_T, World);

    uint64_t num_locreads = lfd->local_count();
    myoffset = 0;

    MPI_Exscan(&num_locreads, &myoffset, 1, MPI_UINT64_T, MPI_SUM, World);

    uint64_t charbuf_totsend; /* will store total number of chars this processor will send */
    uint64_t read_totsend = 0; /* will store total number of reads this processor will send */

    std::vector<int> read_sendcounts(nprocs, 0); /* sendcounts for read indices and read lengths (MPI_Alltoallv::sendcounts) */

    /* In our MPI_Alltoallv call, we will have to store the reads we send (for both the char buffer
     * and associated metadata) packed together by processor. Because the local reads we are sending
     * are assigned more or less randomly, we have to first go through each local read index and
     * determine exactly which processors are receivers of our local reads, so that we can then
     * pack all reads of a particular processor together in the final MPI_Alltoallv send buffers */

    std::vector<std::vector<uint64_t>> i_read_send_idxs(nprocs); /* intra-processor read global indices to send */
    std::vector<std::vector<uint64_t>> i_read_send_lens(nprocs); /* intra-processor read lengths to send */

    std::unordered_map<uint64_t, std::tuple<ushort, uint64_t>> lfd_info_cache; /* cache used read lengths and offsets for later */

    for (uint64_t i = 0; i < num_locreads; ++i) {
        int dest = read2procs_global[i+myoffset];
        /* dest == -1 means the read is not participating in contig creation */
        if (dest >= 0 && myrank != dest) {
            ushort readlen;
            uint64_t start_offset, end_offset;
            lfd->get_sequence(i, readlen, start_offset, end_offset);
            lfd_info_cache[i] = std::make_tuple(readlen, start_offset);
            i_read_send_idxs[dest].push_back(i+myoffset);
            i_read_send_lens[dest].push_back(readlen);
            ++read_sendcounts[dest];
            ++read_totsend;
        }
    }

    std::vector<int> read_sdispls(nprocs, 0); /* read index/length buffer processor send displacements (MPI_Alltoallv::sdispls) */

    std::vector<uint64_t> read_send_idxs(read_totsend); /* read index send buffer (MPI_Alltoallv::sendbuf) */
    std::vector<uint64_t> read_send_lens(read_totsend); /* read length send buffer (MPI_Alltoallv::sendbuf) */

    std::partial_sum(read_sendcounts.begin(), read_sendcounts.end()-1, read_sdispls.begin()+1);

    std::vector<int> charbuf_sendcounts(nprocs, 0); /* sendcounts for packed char buffer (MPI_Alltoallv::sendcounts) */
    std::vector<int> charbuf_sdispls(nprocs, 0); /* send displacements for packed char buffer (MPI_Alltoallv::sdispls) */

    charbuf_totsend = 0;

    /* intra-processor read lengths/indices are ready to be copied into actual send buffers.
     * loop jamming used also to calculate how many chars are being sent to each process */
    for (int i = 0; i < nprocs; ++i) {
        std::copy(i_read_send_idxs[i].begin(), i_read_send_idxs[i].end(), read_send_idxs.begin() + read_sdispls[i]);
        std::copy(i_read_send_lens[i].begin(), i_read_send_lens[i].end(), read_send_lens.begin() + read_sdispls[i]);
        uint64_t rslen = std::accumulate(i_read_send_lens[i].begin(), i_read_send_lens[i].end(), static_cast<uint64_t>(0));
        charbuf_sendcounts[i] = rslen;
        charbuf_totsend += rslen;
    }

    /* TODO :: Can the above loop also include char buffer loading?? */

    std::partial_sum(charbuf_sendcounts.begin(), charbuf_sendcounts.end()-1, charbuf_sdispls.begin()+1);

    std::vector<int> read_recvcounts(nprocs); /* recvcounts for read indices/lengths (MPI_Alltoallv::recvcounts) */
    std::vector<int> read_rdispls(nprocs, 0); /* rdispls for read indices/lengths (MPI_Alltoallv::rdispls) */

    std::vector<int> charbuf_recvcounts(nprocs); /* recvcounts for packed char buffer (MPI_Alltoallv::recvcounts) */
    std::vector<int> charbuf_rdispls(nprocs, 0); /* rdispls for packed char buffer (MPI_Alltoalv::rdispls) */

    /* MPI_Alltoalls so that each processor knows how much it is receiving and from whom */
    MPI_Alltoall(read_sendcounts.data(), 1, MPI_INT, read_recvcounts.data(), 1, MPI_INT, World);
    MPI_Alltoall(charbuf_sendcounts.data(), 1, MPI_INT, charbuf_recvcounts.data(), 1, MPI_INT, World);

    uint64_t read_totrecv = std::accumulate(read_recvcounts.begin(), read_recvcounts.end(), static_cast<uint64_t>(0)); /* total number of reads to receive on this processor */
    std::partial_sum(read_recvcounts.begin(), read_recvcounts.end()-1, read_rdispls.begin()+1);

    uint64_t charbuf_totrecv = std::accumulate(charbuf_recvcounts.begin(), charbuf_recvcounts.end(), static_cast<uint64_t>(0)); /* total number of chars to receive on this processor */
    std::partial_sum(charbuf_recvcounts.begin(), charbuf_recvcounts.end()-1, charbuf_rdispls.begin()+1);

    std::vector<uint64_t> read_recv_idxs(read_totrecv); /* received read indices buffer (MPI_Alltoallv::recvbuf) */
    std::vector<uint64_t> read_recv_lens(read_totrecv); /* received read lengths buffer (MPI_Alltoallv::recvbuf) */

    /* MPI_Alltoallvs so that each processor receives the global read indices that are to be sent
     * here, and the read lengths as well */
    MPI_Alltoallv(read_send_idxs.data(), read_sendcounts.data(), read_sdispls.data(), MPI_UINT64_T, read_recv_idxs.data(), read_recvcounts.data(), read_rdispls.data(), MPI_UINT64_T, World);
    MPI_Alltoallv(read_send_lens.data(), read_sendcounts.data(), read_sdispls.data(), MPI_UINT64_T, read_recv_lens.data(), read_recvcounts.data(), read_rdispls.data(), MPI_UINT64_T, World);

    char *charbuf_send = new char[charbuf_totsend+1](); /* char buffer to send all the reads (MPI_Alltoallv::sendbuf) */
    char *charbuf_send_ptr = charbuf_send;

    const char *lfd_buffer = lfd->buffer();

    ushort len;
    uint64_t offset;

    /* create the packed char buffer to send */
    for (auto itr = read_send_idxs.begin(); itr != read_send_idxs.end(); ++itr) {
        auto cached = lfd_info_cache[*itr - myoffset];
        len = std::get<0>(cached), offset = std::get<1>(cached);
        std::strncpy(charbuf_send_ptr, lfd_buffer + offset, len);
        charbuf_send_ptr += len;
    }

    char *charbuf_recv = new char[charbuf_totrecv+1]();

    MPI_Alltoallv(charbuf_send, charbuf_sendcounts.data(), charbuf_sdispls.data(), MPI_CHAR, charbuf_recv, charbuf_recvcounts.data(), charbuf_rdispls.data(), MPI_CHAR, World);

    delete [] charbuf_send;

    charbuf_read_idxs = read_recv_idxs;
    charbuf_read_lengths = read_recv_lens;

    charbuf_read_offsets.resize(charbuf_read_lengths.size());
    charbuf_read_offsets[0] = 0;

    std::partial_sum(read_recv_lens.begin(), read_recv_lens.end()-1, charbuf_read_offsets.begin()+1);

    *charbuf = charbuf_recv;
    return read_totrecv;
}


#endif
