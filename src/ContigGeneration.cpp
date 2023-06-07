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

    degrees.DebugPrint();

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

std::vector<std::string> GenerateContigs(CT<Overlap>::PSpParMat& S, const DnaBuffer& mydna, DistributedFastaData& dfd)
{
    std::vector<std::string> contigs;

    /*
     * Compute a mapping of read IDs to contig IDs. The function @GetRead2Contigs
     * is responsible for computing the contig sets from the string graph.
     */
    CT<int64_t>::PDistVec assignments(S.getcommgrid());
    int64_t numcontigs = GetRead2Contigs(S, assignments);

    /*
     * Compute contig sizes.
     */
    auto contigsizes = GetContigSizes(assignments, numcontigs, dfd);

    return contigs;
}
