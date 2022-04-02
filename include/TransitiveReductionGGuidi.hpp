#ifndef TRANSITIVE_REDUCTION_H_
#define TRANSITIVE_REDUCTION_H_

#include "../include/ReadOverlap.hpp"

#include <sys/time.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <sstream>

struct InvalidSRing : unary_function<ReadOverlap, ReadOverlap>
{
    bool operator() (const ReadOverlap& x) { return x.is_invalid(); }
};

struct TransposeSRing : unary_function <ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (const ReadOverlap& x) const
    {
        ReadOverlap xT = x;

        xT.b[0] = x.l[1] - x.e[1];
        xT.e[0] = x.l[1] - x.b[1];
        xT.b[1] = x.l[0] - x.e[0];
        xT.e[1] = x.l[0] - x.b[0];

        xT.l[0] = x.l[1];
        xT.l[1] = x.l[0];

        xT.rc = x.rc;
        xT.transpose = !x.transpose;

        xT.sfx  = x.sfxT;
        xT.sfxT = x.sfx;
        xT.dir  = x.dirT;
        xT.dirT = x.dir;

        return xT;
    }
};

struct PlusFuzzSRing : unary_function<ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& x) const
    {
        ReadOverlap fuzzed = x;
        fuzzed.sfx += FUZZ;
        fuzzed.sfxT += FUZZ;

        return fuzzed;
    }
};

struct TransitiveSelection : binary_function<ReadOverlap, OverlapPath, bool>
{
    bool operator() (const ReadOverlap& r, const OverlapPath& n) const
    {
        int dir = r.dir;

        if (dir == -1) return false;

        return (r.sfx >= n.sfx[dir]);
    }
};

struct TransitiveRemoval : binary_function<ReadOverlap, bool, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& x, const bool& y)
    {
        if (!y) x.dir = -1;
        return x;
    }
};

struct ZeroPrune { bool operator() (const bool& x) { return false; } };
struct BoolPrune { bool operator() (const bool& x) { return !x; } };

OverlapPath opmin(const OverlapPath& e1, const OverlapPath& e2)
{
    OverlapPath e = OverlapPath();

    for (int i = 0; i < 4; ++i)
        e.sfx[i] = std::min(e1.sfx[i], e2.sfx[i]);

    return e;
}

struct MinPlusSR
{
    static OverlapPath id() { return OverlapPath(); }
    static bool returnedSAID() { return false; }
    static MPI_Op mpi_op() { return MPI_MIN; }

    static OverlapPath add(const OverlapPath& e1, const OverlapPath& e2)
    {
        return opmin(e1, e2);
    }

    static OverlapPath multiply(const ReadOverlap& e1, const ReadOverlap& e2)
    {
        OverlapPath e = OverlapPath();

        int t1, t2, h1, h2;

        if (!e1.arrows(t1, h1) || !e2.arrows(t2, h2))
            return e;

        if (t2 == h1)
            return e;

        e.sfx[2*t1 + h2] = e1.sfx + e2.sfx;

        return e;
    }

    static void axpy(ReadOverlap a, const ReadOverlap& x, OverlapPath& y)
    {
        y = opmin(y, multiply(a, x));
    }
};

/* TR main function */
#define MAXITER 15

void TransitiveReductionOld(PSpMat<ReadOverlap>::MPI_DCCols& B, TraceUtils tu)
{

    PSpMat<ReadOverlap>::MPI_DCCols RT = B; /* copies everything */
    RT.Transpose();
    RT.Apply(TransposeSRing()); /* flips all the coordinates */

    if (!(RT == R)) /* symmetricize */
    {
        R += RT;
    }    

    R.ParallelWriteMM("overlap-graph-symmetric.mm", true, ReadOverlapExtraHandler());

    /* implicitly will call OverlapPath(const ReadOverlap& e) constructor */
    PSpMat<OverlapPath>::MPI_DCCols P = R; /* P is a copy of R now but it is going to be the "power" matrix to be updated over and over */

    /* create an empty boolean matrix using the same proc grid as R */
    PSpMat<bool>::MPI_DCCols T(R.getcommgrid()); /* T is going to store transitive edges to be removed from R in the end */

    /* create a copy of R and add a FUZZ constant to it so it's more robust to error in the sequences/alignment */ 
    PSpMat<ReadOverlap>::MPI_DCCols F = R;
    F.Apply(PlusFuzzSRing());

    tu.print_str("Matrix R, i.e. AAt post alignment and before transitive reduction: ");
    R.PrintInfo();

    uint nnz, prev;

    uint count = 0;     
    uint countidle = 0; 

    /* Gonna iterate on T until there are no more transitive edges to remove */
    do
    {
        prev = T.getnnz();

        /* Computer N (neighbor matrix)
         * N = P*R
         */
        PSpMat<dibella::CommonKmers>::MPI_DCCols N = Mult_AnXBn_DoubleBuff<MinPlusSR_t, OverlapPath, PSpMat<OverlapPath>::DCCols>(P, R);
        // N.Prune(InvalidSRing(), true); // GGGG this is sketchy to me, let's discuss it

        P = N;

        OverlapPath id;

        /* GGGG: Ii is going to be true is the Ni dir entry corresponding to the Ri dir entry is non-zero, that is
            mark true edges that should be removed from R eventually because transitive.
            I would call this semiring something like "GreaterThenTR", but that's minor.
            */
        PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, SpDCCols<int64_t, bool>>(F, N, TransitiveSelection(), false, id);
        
        /* GGGG: have no idea why this stuff is happening here
            to make sure every transitive edge is correctly removed in the upper/lower triangular entry, too
            this would happen naturally on an error-free dataset 
            */

        /* I has some 1 entries, now IT has them, too */
        PSpMat<bool>::MPI_DCCols IT = I;
        /* IT is not transpose, so IT has some 1 entries != from the 1 entries of I */
        IT.Transpose();

        /* GGGG: if you element-wise multiply an entry that is set to an entry in I times an entry set to 1 in IT but 0 in I, isn't the entry in I staying 0 rather becoming 1? 
            Don't you want to do somethine like:
        
            bool isLogicalNot = false;
            bool bId = false;
            I = EWiseApply<bool, SpDCCols<int64_t, bool>>(I, IT, bind2nd, isLogicalNot, bId); 

            where bind2nd would store 1 in I[i, j] if IT[i, j] is 1 and do nothing otherwise (i.e., if I[i, j] is 1, it's gonna stay that way even if IT[i, j] is 0)
            */
        I.EWiseMult(IT, false); // GGGG: this operation doesn't make sense to me, missing something? (see comment right above)

        bool isLogicalNot = true;
        bool bId = false;

        /* GGGG: this semiring/lambda does NAND, right? But shouldn't this only return 0 if both are 0? I'm not sure I understand why we are returning 0 where they are both 1? 
           If T(i,j) = 1 and I(i,j) = 1 why the resulting T(i,j) should be 0? Don't we wanna remove that?
           this should be an OR, not a NAND: return (x || y);
        */
        T = EWiseApply<bool, SpDCCols<int64_t, bool>>(T, I, [](bool x, bool y) { return (x || y); }, isLogicalNot, bId); 
        cur = T.getnnz();

        /* GGGG:
            (prev - cur) < 100 means that there are less than 100 nonzeros of difference between T at timestep i and T at timestep i-1
            check_phase_counter is inizialied to 0
            if there are less than 100 differences, then increment check_phase_counter by 1
            once we increment it the first time, then it keeps incrementing until it reaches 5 and the do-while loop is gonna stop,
            if we don't trigger the (prev - cur) < 100 condition, the do-while loop will stop regardless after 10 iterations.     
            */

        /* GGGG: this is probably unnecessarily complicated, let's do something like 'if nonzeros in T don't change for X iterations, then stop' */
        // if (check_phase_counter > 0 || (prev - cur) < 100)
        // {
        //     check_phase_counter++;
        // }

        // full_counter++;

        /* GGGG: if nonzeros in T don't change for MAXITER iterations, then exit the loop 
            If in the test run this creates isses, e.g. it doesn't terminate, let's put a stricter exit condition */
        if(cur == prev)
        {
            countidle++; 
        }
        else if(countidle > 0)
        {
            countidle = 0; // GGGG: reset the counter if cur != prev in this iteration but countidle > 0
        }

        count++; // GGGG: this just keeps track of the total number of iterations but doesn't do anything about the termination condition

    } while (countidle < MAXITER)
    // } while (check_phase_counter < 5 && full_counter < 10); 

    std::stringstream iss;
    iss << "the transitive reduction did " << count << " iterations\n";
    tu.print_str(iss.str());

    bool isLogicalNot = false;
    bool bId = false;

    R = EWiseApply<ReadOverlap, SpDCCols<int64_t, ReadOverlap>>(R, T, TransitiveRemoval(), isLogicalNot, bId);
    R.Prune(InvalidSRing());

    R.ParallelWriteMM("string-graph.mm", true, ReadOverlapExtraHandler());

    tu.print_str("Matrix S, i.e. AAt post transitive reduction: ");
    R.PrintInfo();
}

#endif
