#include "ContigGeneration.hpp"
#include "CC.hpp"

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

CT<int64_t>::PDistVec GetContigSizes(const CT<int64_t>::PDistVec& assignments, int64_t numcontigs, DistributedFastaData& dfd)
{
    auto commgrid = assignments.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    std::vector<int64_t> local_component_sizes(numcontigs, 0);
    std::vector<int64_t> local_components = assignments.GetLocVec();

    int64_t local_size = local_components.size();

    for (int64_t i = 0; i < local_size; ++i)
    {
        ++local_component_sizes[local_components[i]];
    }

    int avgsize = numcontigs / nprocs;
    int lastsize = numcontigs - (avgsize * (nprocs - 1));

    std::vector<int> recvcounts(nprocs, avgsize);
    recvcounts.back() = lastsize;

    int mysize = (myrank != nprocs-1)? avgsize : lastsize;

    std::vector<int64_t> fill_vec_components(mysize);

    MPI_Reduce_scatter(local_component_sizes.data(), fill_vec_components.data(), recvcounts.data(), MPI_INT64_T, MPI_SUM, comm);

    return CT<int64_t>::PDistVec(fill_vec_components, commgrid);
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

    assignments.DebugPrint();

    /*
     * Compute contig sizes.
     */
    CT<int64_t>::PDistVec contigsizes = GetContigSizes(assignments, numcontigs, dfd);

    contigsizes.DebugPrint();

    return contigs;
}
