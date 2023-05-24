#include "TransitiveReduction.hpp"

std::unique_ptr<CT<Overlap>::PSpParMat> TransitiveReduction(CT<Overlap>::PSpParMat R)
{
    auto commgrid = R.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    /*
     * @R is upper triangular, but we need it to be symmetric to compute the
     * transitive reduction.
     */

    CT<Overlap>::PSpParMat RT = R;
    RT.Transpose();
    RT.Apply(Overlap::Transpose()); /* This flips all query and target info encoded in the nonzeros */

    R += RT; /* symmetricize */

    CT<Overlap>::PSpParMat P = R; /* P is a copy of R now but it is going to be the "power" matrix to be updated over and over */

    uint64_t nrow = R.getnrow();
    uint64_t ncol = R.getncol();

    CT<uint64_t>::PDistVec initvec(commgrid, nrow, static_cast<uint64_t>(0));
    CT<int>::PSpParMat T(nrow, ncol, initvec, initvec, static_cast<int>(0), false);

    /*
     * Create a copy of R and add a FUZZ constant to it so it's more robust to error in the sequences/alignment.
     */
    CT<Overlap>::PSpParMat F = R;
    F.Apply(PlusFuzzSRing());

    uint64_t cur, prev;
    uint64_t count = 0;
    uint64_t countidle = 0;

    bool isLogicalNot = false;
    bool bId = false;

    /*
     * Iterate on T until there are no more transitive edges to remove.
     */
    do
    {
        prev = T.getnnz();
        CT<Overlap>::PSpParMat N = Mult_AnXBn_DoubleBuff<MinPlusSR, Overlap, CT<Overlap>::PSpDCCols>(P, R);

        N.Prune(NoPathSRing(), true);

        P = N;
        Overlap id;

        /* GGGG:
         * Ii is going to be true is the Ni dir entry corresponding to the Ri dir entry is non-zero, that is
         * mark true edges that should be removed from R eventually because transitive.
         */
        CT<bool>::PSpParMat I = EWiseApply<bool, CT<bool>::PSpDCCols>(F, N, GreaterThanSR(), false, id);

        I.Prune(BoolPrune(), true); /* GGGG: this is needed to remove entries in N that were smaller than F and thus resulted in an actual 0 in F */

        /* GGGG:
         * make sure every transitive edge is correctly removed in the upper/lower triangular entry, too
         * this would happen naturally on an error-free dataset.
         */

        /* I has some 1 entries, now IT has them, too */
        CT<bool>::PSpParMat IT = I;
        IT.Transpose();

        if (!(IT == I)) I += IT; /* symmetricize */

        T += I; /* GGGG: add new transitive edges to T */

        cur = T.getnnz();

        count++; /* GGGG: this just keeps track of the total number of iterations but doesn't do anything about the termination condition */

    }  while (cur != prev);

    isLogicalNot = true; /* GGGG: I want the ones to be removed to be set to false, so I use the logical negation */
    bId = true; /* GGGG: this is critical, if this would be false, everything that would survive the EWIseApply would have dir == -1 as well and it's wrong
                    A non-transitive edge should keep its original direction! */
    R = EWiseApply<Overlap, CT<Overlap>::PSpDCCols>(R, T, TransitiveRemoval(), isLogicalNot, 1);

    R.Prune(InvalidSRing(), true);

    return std::make_unique<CT<Overlap>::PSpParMat>(R);
}

Overlap opmin(const Overlap& e1, const Overlap& e2)
{
    Overlap e = Overlap();

    for (int i = 0; i < 4; ++i)
        e.suffix_paths[i] = std::min(e1.suffix_paths[i], e2.suffix_paths[i]);

    return e;
}
