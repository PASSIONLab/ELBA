#include "ContigGeneration.hpp"
#include "CC.hpp"
#include "Logger.hpp"
#include <unordered_set>

/*
 * @func GetRead2Contigs     determines which reads belong to which contigs.
 *
 * @param S                  distributed string graph.
 * @param assignments        computed assignments vector.
 *
 * @description
 * Contig branching points are by CombBLAS row reduction (to find vertex degrees), and then
 * selecting rows more than 2 nonzeros. These are temporarily zeroed from matrix so that
 * connected components will find which reads belong to which contigs.
 *
 * @return the number of contigs found.
 */
int64_t GetRead2Contigs(CT<Overlap>::PSpParMat& S, CT<int64_t>::PDistVec& assignments)
{
    /*
     * Create a copy of string graph that mimics that nonzero pattern
     * of the graph, but stores int64_t instead of Overlap as nonzero.
     */
    CT<int64_t>::PSpParMat A = S;
    CT<bool>::PSpParMat D = A;

    /*
     * Row reduction using addition gives us a distributed vector
     * of vertex degrees, because the input matrix is symmetric.
     */
    CT<int64_t>::PDistVec degrees(A.getcommgrid());
    D.Reduce(degrees, Row, std::plus<int64_t>(), static_cast<int64_t>(0));

    /*
     * Find the indices of the degrees vector that have entries greater than 2.
     * Those correspond to the branching vertices.
     */
    CT<int64_t>::PDistVec branches = degrees.FindInds([](const int64_t& v) { return v > 2; });
    int64_t numbranches = branches.TotalLength();

    /*
     * Remove the branching vertices from the shadow copy of the string matrix.
     */
    A.PruneFull(branches, branches);

    /*
     * Compute connected components on the branchless string matrix.
     */
    int64_t numcontigs;
    assignments = CC(A, numcontigs);

    return numcontigs;
}

std::vector<std::tuple<int64_t, int64_t>> GetContigSizes(const CT<int64_t>::PDistVec& assignments, int64_t numcontigs, DistributedFastaData& dfd)
{
    auto commgrid = assignments.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    /*
     * @local_component_sizes is a local vector measuring the number of reads (stored locally)
     * on each contig. Its range is across the entire global set of contigs, but only counts
     * the contributors from my processor.
     */
    std::vector<int64_t> local_component_sizes(numcontigs, 0);
    std::vector<int64_t> local_components = assignments.GetLocVec();

    int64_t local_size = local_components.size();

    /*
     * @local_components[i] is the contig ID of the ith local read. Find
     * its bucket and increment its counter (which is initialized to zero).
     */
    for (int64_t i = 0; i < local_size; ++i)
    {
        ++local_component_sizes[local_components[i]];
    }

    /*
     * Compute and distribute the global contig sizes (number of reads per contig) using
     * a reduce-scatter collective. The results are distributed linearly to each processor.
     * The local partition of contig sizes is stored in @locsizes.
     */

    int avgsize = numcontigs / nprocs;
    int lastsize = numcontigs - (avgsize * (nprocs - 1));

    std::vector<int> recvcounts(nprocs, avgsize);
    recvcounts.back() = lastsize;

    int64_t mysize = (myrank != nprocs-1)? avgsize : lastsize;

    std::vector<int64_t> locsizes(mysize);

    MPI_Reduce_scatter(local_component_sizes.data(), locsizes.data(), recvcounts.data(), MPI_INT64_T, MPI_SUM, comm);

    int64_t lengthuntil;
    MPI_Exscan(&mysize, &lengthuntil, 1, MPI_INT64_T, MPI_SUM, comm);
    if (myrank == 0) lengthuntil = 0;

    /*
     * Construct vector of tuples whose entries are (global contig ID, contig size),
     * where we only store contigs that have at least 2 reads.
     */
    std::vector<std::tuple<int64_t, int64_t>> sendbuf;
    for (int64_t i = 0; i < mysize; ++i)
        if (locsizes[i] >= 2)
            sendbuf.emplace_back(i + lengthuntil, locsizes[i]);

    std::vector<int> displs(nprocs);

    std::fill(recvcounts.begin(), recvcounts.end(), static_cast<int64_t>(0));
    recvcounts[myrank] = sendbuf.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);
    std::exclusive_scan(recvcounts.begin(), recvcounts.end(), displs.begin(), static_cast<int64_t>(0));

    int64_t num_used_contigs = recvcounts.back() + displs.back();

    std::vector<std::tuple<int64_t, int64_t>> result(num_used_contigs);

    MPI_Allgatherv(sendbuf.data(), recvcounts[myrank], MPIType<std::tuple<int64_t, int64_t>>(), result.data(), recvcounts.data(), displs.data(), MPIType<std::tuple<int64_t, int64_t>>(), comm);

    std::sort(result.begin(), result.end(), [](const auto& a, const auto& b) { return std::get<1>(a) > std::get<1>(b); });

    return result;
}

std::vector<int64_t> ImposeMyReadDistribution(CT<int64_t>::PDistVec& assignments, DistributedFastaData& dfd)
{
    auto index = dfd.getindex();
    auto commgrid = index.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    int64_t orig_offset = assignments.LengthUntil();
    int64_t orig_size = assignments.LocArrSize();

    std::vector<int64_t> orig_vector = assignments.GetLocVec();

    std::vector<int> sendcounts(nprocs, 0);
    std::vector<int> recvcounts(nprocs);

    for (int64_t i = 0; i < orig_size; ++i)
    {
        int owner = index.getreadowner(i+orig_offset);
        ++sendcounts[owner];
    }

    MPI_Alltoall(sendcounts.data(), 1, MPI_INT, recvcounts.data(), 1, MPI_INT, comm);

    std::vector<int> sdispls(nprocs);
    std::vector<int> rdispls(nprocs);

    std::exclusive_scan(sendcounts.begin(), sendcounts.end(), sdispls.begin(), 0);
    std::exclusive_scan(recvcounts.begin(), recvcounts.end(), rdispls.begin(), 0);

    int64_t totrecv = recvcounts.back() + rdispls.back();

    std::vector<int64_t> new_vector(totrecv);

    MPI_Alltoallv(orig_vector.data(), sendcounts.data(), sdispls.data(), MPI_INT64_T, new_vector.data(), recvcounts.data(), rdispls.data(), MPI_INT64_T, comm);

    return new_vector;
}

std::vector<int64_t> GetLocalProcAssignments(CT<int64_t>::PDistVec& assignments, std::vector<std::tuple<int64_t, int64_t>>& contigsizes, DistributedFastaData& dfd)
{
    auto index = dfd.getindex();
    auto commgrid = index.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    int64_t num_used_contigs = contigsizes.size();
    std::vector<int64_t> small_proc_assignments(num_used_contigs);
    std::vector<int64_t> small_large_map(num_used_contigs);

    if (myrank == 0)
    {
        std::vector<int64_t> sums(nprocs, 0);
        std::vector<std::vector<int64_t>> partitions(nprocs);

        for (int64_t i = 0; i < num_used_contigs; ++i)
        {
            int where = std::distance(sums.begin(), std::min_element(sums.begin(), sums.end()));
            sums[where] += std::get<1>(contigsizes[i]);
            small_large_map[i] = std::get<0>(contigsizes[i]);
            partitions[where].push_back(i);
        }

        for (int i = 0; i < nprocs; ++i)
            for (auto itr = partitions[i].begin(); itr != partitions[i].end(); ++itr)
                small_proc_assignments[*itr] = i;
    }

    MPI_Bcast(small_proc_assignments.data(), static_cast<int>(num_used_contigs), MPI_INT64_T, 0, comm);
    MPI_Bcast(small_large_map.data(), static_cast<int>(num_used_contigs), MPI_INT64_T, 0, comm);

    std::unordered_map<int64_t, int64_t> large_small_map;

    for (auto itr = small_large_map.begin(); itr != small_large_map.end(); ++itr)
        large_small_map[*itr] = (itr - small_large_map.begin());

    int64_t lengthuntil = index.getmyreaddispl();
    int64_t nlocreads = index.getmyreadcount();

    std::vector<int64_t> local_assignments = ImposeMyReadDistribution(assignments, dfd);
    std::vector<int64_t> local_proc_assignments(nlocreads, -1);

    for (auto itr = local_assignments.begin(); itr != local_assignments.end(); ++itr)
        if (large_small_map.find(*itr) != large_small_map.end())
            local_proc_assignments[itr - local_assignments.begin()] = small_proc_assignments[large_small_map[*itr]];

    return local_proc_assignments;
}

std::unordered_map<int64_t, std::string> GetInducedReadSequences(const std::vector<int64_t> local_proc_assignments, const DnaBuffer& mydna, DistributedFastaData& dfd)
{
    auto index = dfd.getindex();
    auto commgrid = index.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    int64_t myreaddispl = index.getmyreaddispl();

    std::vector<MPI_Count_type> item_sendcnts(nprocs, 0);
    std::vector<MPI_Count_type> item_recvcnts(nprocs);
    std::vector<MPI_Displ_type> item_sdispls(nprocs);
    std::vector<MPI_Displ_type> item_rdispls(nprocs);

    std::vector<MPI_Count_type> char_sendcnts(nprocs, 0);
    std::vector<MPI_Count_type> char_recvcnts(nprocs);
    std::vector<MPI_Displ_type> char_sdispls(nprocs);
    std::vector<MPI_Displ_type> char_rdispls(nprocs);

    std::vector<int64_t> globalid_sendbuf;
    std::vector<int64_t> globalid_recvbuf;
    std::vector<std::vector<int64_t>> globalid_buckets(nprocs);
    int64_t item_totsend, item_totrecv;

    std::vector<int64_t> readlen_sendbuf;
    std::vector<int64_t> readlen_recvbuf;

    std::vector<char> char_sendbuf;
    std::vector<char> char_recvbuf;
    int64_t char_totsend, char_totrecv;

    item_totsend = char_totsend = 0;

    for (int64_t i = 0; i < local_proc_assignments.size(); ++i)
    {
        int dest = local_proc_assignments[i];

        if (dest >= 0)
        {
            item_totsend++;
            item_sendcnts[dest]++;

            char_totsend += mydna[i].size();
            char_sendcnts[dest] += mydna[i].size();

            globalid_buckets[dest].push_back(i + myreaddispl);
        }
    }

    MPI_ALLTOALL(char_sendcnts.data(), 1, MPI_COUNT_TYPE, char_recvcnts.data(), 1, MPI_COUNT_TYPE, comm);

    std::exclusive_scan(char_sendcnts.begin(), char_sendcnts.end(), char_sdispls.begin(), static_cast<MPI_Displ_type>(0));
    std::exclusive_scan(char_recvcnts.begin(), char_recvcnts.end(), char_rdispls.begin(), static_cast<MPI_Displ_type>(0));

    assert(char_totsend == char_sendcnts.back() + char_sdispls.back());
    char_totrecv = char_recvcnts.back() + char_rdispls.back();

    char_sendbuf.reserve(char_totsend);
    char_recvbuf.resize(char_totrecv);

    MPI_ALLTOALL(item_sendcnts.data(), 1, MPI_COUNT_TYPE, item_recvcnts.data(), 1, MPI_COUNT_TYPE, comm);

    std::exclusive_scan(item_sendcnts.begin(), item_sendcnts.end(), item_sdispls.begin(), static_cast<MPI_Displ_type>(0));
    std::exclusive_scan(item_recvcnts.begin(), item_recvcnts.end(), item_rdispls.begin(), static_cast<MPI_Displ_type>(0));

    assert(item_totsend == item_sendcnts.back() + item_sdispls.back());
    item_totrecv = item_recvcnts.back() + item_rdispls.back();

    readlen_sendbuf.reserve(item_totsend);
    readlen_recvbuf.resize(item_totrecv);

    globalid_sendbuf.reserve(item_totsend);
    globalid_recvbuf.resize(item_totrecv);

    for (auto itr = globalid_buckets.begin(); itr != globalid_buckets.end(); ++itr)
    {
        globalid_sendbuf.insert(globalid_sendbuf.end(), itr->begin(), itr->end());

        for (size_t i = 0; i < itr->size(); ++i)
        {
            int64_t localid = (*itr)[i] - myreaddispl;
            auto seq = mydna[localid].ascii();

            char_sendbuf.insert(char_sendbuf.end(), seq.begin(), seq.end());
            readlen_sendbuf.push_back(seq.size());
        }
    }

    MPI_ALLTOALLV(globalid_sendbuf.data(), item_sendcnts.data(), item_sdispls.data(), MPI_INT64_T,
                  globalid_recvbuf.data(), item_recvcnts.data(), item_rdispls.data(), MPI_INT64_T, comm);

    MPI_ALLTOALLV(readlen_sendbuf.data(), item_sendcnts.data(), item_sdispls.data(), MPI_INT64_T,
                  readlen_recvbuf.data(), item_recvcnts.data(), item_rdispls.data(), MPI_INT64_T, comm);

    MPI_ALLTOALLV(char_sendbuf.data(), char_sendcnts.data(), char_sdispls.data(), MPI_CHAR,
                  char_recvbuf.data(), char_recvcnts.data(), char_rdispls.data(), MPI_CHAR, comm);

    std::unordered_map<int64_t, std::string> read_lookup_table;

    auto itr = char_recvbuf.begin();

    for (int64_t i = 0; i < item_totrecv; ++i)
    {
        int64_t globalid = globalid_recvbuf[i];
        size_t readlen = readlen_recvbuf[i];
        std::string seq(itr, itr + readlen);
        read_lookup_table[globalid] = seq;
        itr += readlen;
    }

    return read_lookup_table;
}

char comp(char c)
{
    switch (c)
    {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
    }
    return '\0';
}

std::vector<std::string> GenerateContigs(CT<Overlap>::PSpParMat& S, const DnaBuffer& mydna, DistributedFastaData& dfd)
{
    auto index = dfd.getindex();
    auto commgrid = index.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();
    std::vector<std::string> contigs;

    /*
     * Compute a mapping of read IDs to contig IDs. The function @GetRead2Contigs
     * is responsible for computing the contig sets from the string graph.
     */
    CT<int64_t>::PDistVec assignments(commgrid);
    int64_t numcontigs = GetRead2Contigs(S, assignments);

    /*
     * Compute contig sizes.
     */
    auto contigsizes = GetContigSizes(assignments, numcontigs, dfd);

    /*
     * Compute local read ID to processor assignments.
     */
    std::vector<int64_t> local_proc_assignments = GetLocalProcAssignments(assignments, contigsizes, dfd);

    CT<int64_t>::PDistVec proc_assignments(local_proc_assignments, commgrid);

    std::vector<int64_t> local_contig_read_ids; /* mapping of local graph indices (`contig_chains`) to original global read indices. */
    CT<Overlap>::PSpDCCols contig_chains_derived = S.InducedSubgraphs2Procs(proc_assignments, local_contig_read_ids);
    CT<Overlap>::PSpCCols contig_chains(contig_chains_derived);
    contig_chains.Transpose();

    auto read_lookup_table = GetInducedReadSequences(local_proc_assignments, mydna, dfd);

    int64_t myreaddispl = index.getmyreaddispl();
    int64_t numreads = contig_chains.getnrow();

    Logger logger(commgrid);

    assert(numreads == contig_chains.getncol());
    assert(numreads == read_lookup_table.size());
    assert(numreads == local_contig_read_ids.size());

    auto csc = contig_chains.GetCSC();

    if (!csc) return contigs;

    assert(numreads >= 2);

    std::vector<bool> visited(numreads, false);
    std::unordered_set<int64_t> used_roots;

    int64_t contig_id = 0;

    for (int64_t v = 0; v < csc->n; ++v)
    {
        if (csc->jc[v+1] - csc->jc[v] != 1 || used_roots.find(v) != used_roots.end())
            continue;

        std::vector<std::tuple<int64_t, int64_t, bool>> chain;
        int lastdir;

        int64_t cur = v;
        int64_t next, end;

        while (true)
        {
            visited[cur] = true;
            next = csc->jc[cur];
            end = csc->jc[cur + 1];

            while (next < end && visited[csc->ir[next]])
                ++next;

            if (next >= end)
                break;

            Overlap& o = csc->num[next];
            int strand = (o.direction >> 1) & 1;
            chain.emplace_back(local_contig_read_ids[cur], o.suffixT, static_cast<bool>(strand));
            lastdir = o.direction;

            cur = csc->ir[next];
        }

        int64_t readlen = read_lookup_table[local_contig_read_ids[cur]].size();
        chain.emplace_back(local_contig_read_ids[cur], readlen, static_cast<bool>(1 - (lastdir & 1)));

        std::string contig = "";

        for (auto itr = chain.begin(); itr != chain.end(); ++itr)
        {
            int64_t readid = std::get<0>(*itr);
            int64_t prefix = std::get<1>(*itr);
            bool strand = std::get<2>(*itr);

            auto s = read_lookup_table[readid];

            if (strand)
            {
                std::transform(s.cbegin(), s.cend(), s.begin(), comp);
                std::reverse(s.begin(), s.end());
            }

            contig += std::string(s.begin(), s.begin() + prefix);
        }

        contigs.push_back(contig);

        used_roots.insert(cur);
    }

    return contigs;
}
