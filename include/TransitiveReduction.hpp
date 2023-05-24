#ifndef TRANSITIVE_REDUCTION_H_
#define TRANSITIVE_REDUCTION_H_

#include "Overlap.hpp"
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <sstream>

#define FUZZ (1000)

//#define SETITER 2
//#define PDEBUG

/* TR main function */
// #define BECONSERVATIVE
#define TIMINGTR
#ifdef  BECONSERVATIVE
#define MAXITER 5
#endif

struct InvalidSRing : unary_function<Overlap, Overlap>
{
    bool operator() (const Overlap& x) { return (x.direction == -1); }
};

struct NoPathSRing : unary_function<Overlap, Overlap>
{
    bool operator() (const Overlap& x)
    {
        for (int i = 0; i < 4; ++i)
            if (x.suffix_paths[i] < std::numeric_limits<int>::max())
                return false;

        return true;
    }
};

struct PlusFuzzSRing : unary_function<Overlap, Overlap>
{
    Overlap operator() (Overlap& x) const
    {
        Overlap fuzzed = x;
        fuzzed.suffix  += FUZZ;
        fuzzed.suffixT += FUZZ;

        return fuzzed;
    }
};

struct GreaterThanSR : binary_function<Overlap, Overlap, bool>
{
    bool operator() (const Overlap& r, const Overlap& n) const
    {
        int dir = r.direction;

        if (dir == -1) return false;

        return (r.suffix >= n.suffix_paths[dir]);
    }
};

struct TransitiveRemoval : binary_function<Overlap, int, Overlap>
{
    Overlap operator() (Overlap& x, const int& y)
    {
        if (!y)
        {
            x.direction = -1;
        }

        return x;
    }
};

struct ZeroPrune  { bool operator() (const bool& x) { /* If something is a nonzero inherited from R, just prune it! */ return true; } };
struct BoolPrune  { bool operator() (const bool& x) { return !x; } };

Overlap opmin(const Overlap& e1, const Overlap& e2)
{
    Overlap e = Overlap();

    for (int i = 0; i < 4; ++i)
        e.suffix_paths[i] = std::min(e1.suffix_paths[i], e2.suffix_paths[i]);

    return e;
}

struct MinPlusSR
{
    static Overlap id() { return Overlap(); }
    static bool returnedSAID() { return false; }
    static MPI_Op mpi_op() { return MPI_MIN; }

    static Overlap add(const Overlap& e1, const Overlap& e2)
    {
        return opmin(e1, e2);
    }

    static Overlap multiply(const Overlap& e1, const Overlap& e2)
    {
        Overlap e = Overlap();

        int t1, t2, h1, h2;

        if (!e1.arrows(t1, h1) || !e2.arrows(t2, h2))
            return e;

        if (t2 == h1)
            return e;

        e.suffix_paths[2*t1 + h2] = e1.suffix + e2.suffix;

        return e;
    }

    static void axpy(Overlap a, const Overlap& x, Overlap& y)
    {
        y = opmin(y, multiply(a, x));
    }
};

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

//void TransitiveReduction(CT<Overlap>::PSpParMat& R)
//{
//    auto commgrid = R.getcommgrid();
//    int myrank = commgrid->GetRank();
//    int nprocs = commgrid->GetSize();
//    MPI_Comm comm = commgrid->GetWorld();
//
//
//    CT<Overlap>::PSpParMat RT = R; /* copies everything */
//    RT.Transpose();
//    RT.Apply(Overlap::Transpose()); /* swaps query and target coordinates for each nonzero */
//
//    R += RT; /* symmetricize */
//
//    CT<Overlap>::PSpParMat P = R; /* P is a copy of R now but it is going to be the "power" matrix to be updated over and over */
//
//    uint64_t nrow = R.getnrow();
//    uint64_t ncol = R.getncol();
//
//    assert((nrow == ncol));
//
//    CT<uint64_t>::PDistVec initvec(commgrid, nrow, static_cast<uint64_t>(0));
//    CT<int>::PSpParMat T(nrow, ncol, initvec, initvec, static_cast<int>(0), false);
//
//    /* create a copy of R and add a FUZZ constant to it so it's more robust to error in the sequences/alignment */
//    CT<Overlap>::PSpParMat F = R;
//    F.Apply(PlusFuzzSRing());
//
//    tu.print_str("Matrix R, i.e. AAt post alignment and before transitive reduction: ");
//    R.PrintInfo();
//
//    uint64_t cur, prev;
//
//    uint64_t count = 0;
//    uint64_t countidle = 0;
//
//    bool isLogicalNot = false;
//    bool bId = false;
//
//    /* Gonna iterate on T until there are no more transitive edges to remove */
//    do
//    {
//        prev = T.getnnz();
//        CT<Overlap>::PSpParMat N = Mult_AnXBn_DoubleBuff<MinPlusSR, Overlap, CT<Overlap>::PSpDCCols>(P, R);
//
//        N.Prune(NoPathSRing(), true); /* GRGR changed */ /* GGGG: let's discuss this */
//
//        P = N;
//        Overlap id;
//
//        /* GGGG: Ii is going to be true is the Ni dir entry corresponding to the Ri dir entry is non-zero, that is
//            mark true edges that should be removed from R eventually because transitive.
//            I would call this semiring something like "GreaterThenTR", but that's minor.
//            */
//
//        CT<bool>::PSpParMat I = EWiseApply<bool, CT<bool>::PSpDCCols>(F, N, TransitiveSelection(), false, id);
//
//        I.Prune(BoolPrune(), true); /* GGGG: this is needed to remove entries in N that were smaller than F and thus resulted in an actual 0 in F */
//
//        /* GGGG: make sure every transitive edge is correctly removed in the upper/lower triangular entry, too
//            this would happen naturally on an error-free dataset
//            */
//
//        /* I has some 1 entries, now IT has them, too */
//        CT<bool>::PSpParMat IT = I;
//        IT.Transpose();
//
//        if (!(IT == I)) /* symmetricize */
//        {
//            I += IT;
//        }
//
//        T += I; /* GGGG: add new transitive edges to T */
//
//        cur = T.getnnz();
//
//        /* GGGG: if nonzeros in T don't change for MAXITER iterations, then exit the loop
//            If in the test run this creates isses, e.g. it doesn't terminate, let's put a stricter exit condition */
//    #ifdef BECONSERVATIVE
//        if(cur == prev)
//        {
//            countidle++;
//        }
//        else if(countidle > 0)
//        {
//            countidle = 0; // GGGG: reset the counter if cur != prev in this iteration but countidle > 0
//        }
//    #endif
//
//        count++; // GGGG: this just keeps track of the total number of iterations but doesn't do anything about the termination condition
//
//#ifdef SETITER
//    } while (count < SETITER);
//#elif defined(BECONSERVATIVE)
//    } while (countidle < MAXITER);
//#else
//    }  while (cur != prev);
//#endif
//
//    /* GGGG: this is not working as expected! there was a problem in the semiring but it's not fully fixed */
//    R.PrintInfo();
//
//    isLogicalNot = true; /* GGGG: I want the ones to be removed to be set to false, so I use the logical negation */
//    bId = true; /* GGGG: this is critical, if this would be false, everything that would survive the EWIseApply would have dir == -1 as well and it's wrong
//                    A non-transitive edge should keep its original direction! */
//    // TODO R = EWiseApply<ReadOverlap, SpDCCols<int64_t, ReadOverlap>>(R, T, TransitiveRemoval(), isLogicalNot, static_cast<int>(1));
//
//    R.Prune(InvalidSRing(), true);
//}

#endif
