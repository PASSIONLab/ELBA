#include "ContigGeneration.hpp"
#include "CC.hpp"
#include "Logger.hpp"

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

    Logger logger(commgrid);
    for (auto itr = small_proc_assignments.begin(); itr != small_proc_assignments.end(); ++itr)
    {
        logger() << *itr << " ";
    }
    logger.Flush("small_proc_assignments");

    for (auto itr = small_large_map.begin(); itr != small_large_map.end(); ++itr)
    {
        logger() << *itr << " ";
    }
    logger.Flush("small_large_map");

    std::unordered_map<int64_t, int64_t> large_small_map;

    for (auto itr = small_large_map.begin(); itr != small_large_map.end(); ++itr)
        large_small_map[*itr] = (itr - small_large_map.begin());

    for (auto itr = large_small_map.begin(); itr != large_small_map.end(); ++itr)
    {
        logger() << itr->first << ":" << itr->second << ", ";
    }
    logger.Flush("large_small_map");

    int64_t lengthuntil = index.getmyreaddispl();
    int64_t nlocreads = index.getmyreadcount();

    std::vector<int64_t> local_assignments = ImposeMyReadDistribution(assignments, dfd);
    std::vector<int64_t> local_proc_assignments(nlocreads, -1);

    for (auto itr = local_assignments.begin(); itr != local_assignments.end(); ++itr)
    {
        logger() << *itr << " ";
    }
    logger.Flush("local_assignments");

    for (auto itr = local_assignments.begin(); itr != local_assignments.end(); ++itr)
        if (large_small_map.find(*itr) != large_small_map.end())
            local_proc_assignments[itr - local_assignments.begin()] = small_proc_assignments[large_small_map[*itr]];

    return local_proc_assignments;
}

std::vector<std::string> GenerateContigs(CT<Overlap>::PSpParMat& S, const DnaBuffer& mydna, DistributedFastaData& dfd)
{
    auto commgrid = S.getcommgrid();
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

    std::vector<int64_t> local_contig_read_ids;

    CT<Overlap>::PSpDCCols contig_chains_derived = S.InducedSubgraphs2Procs(proc_assignments, local_contig_read_ids);

    CT<Overlap>::PSpCCols contig_chains(contig_chains_derived);
    contig_chains.Transpose();

    return contigs;
}
