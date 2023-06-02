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

std::vector<std::string> GenerateContigs(CT<Overlap>::PSpParMat& S, const DnaBuffer& mydna, DistributedFastaData& dfd)
{
    std::vector<std::string> contigs;

    CT<int64_t>::PDistVec assignments(S.getcommgrid());

    GetRead2Contigs(S, assignments);
    return contigs;
}
