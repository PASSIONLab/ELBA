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

/* TR main function */
// #define BECONSERVATIVE
#define TIMINGTR
#ifdef  BECONSERVATIVE
#define MAXITER 5
#endif

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
        fuzzed.sfx  += FUZZ;
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
        if (!y) 
        {
            x.dir = -1; /* GGGG: This used to be !y and is wrong, we want to removed stuff from R that is set to true in T, not to false */
        }

        return x;
    }
};

struct ZeroPrune  { bool operator() (const bool& x) { /* If something is a nonzero inherited from R, just prune it! */ return true; } }; 
struct BoolPrune  { bool operator() (const bool& x) { return !x; } };

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

void TransitiveReduction(PSpMat<ReadOverlap>::MPI_DCCols& R, TraceUtils tu)
{

    PSpMat<ReadOverlap>::MPI_DCCols RT = R; /* copies everything */
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
    PSpMat<bool>::MPI_DCCols T = R; //(R.getcommgrid()); /* T is going to store transitive edges to be removed from R in the end */
    T.Prune(ZeroPrune()); /* GGGG: there's a better way, TODO for myself */    
    
    /* create a copy of R and add a FUZZ constant to it so it's more robust to error in the sequences/alignment */ 
    PSpMat<ReadOverlap>::MPI_DCCols F = R;
    F.Apply(PlusFuzzSRing());

    tu.print_str("Matrix R, i.e. AAt post alignment and before transitive reduction: ");
    R.PrintInfo();
    tu.print_str("\n");

    uint cur, prev;

    uint count = 0;     
    uint countidle = 0; 
    
    bool isLogicalNot = false;
    bool bId = false;

    double timePR = 0, timeI = 0, timeT = 0, timeTR = 0;

    /* Gonna iterate on T until there are no more transitive edges to remove */
    do
    {
        prev = T.getnnz();      
    #ifdef PDEBUG
        std::cout << prev << " prev" << std::endl; 
    #endif
        /* Computer N (neighbor matrix)
         * N = P*R
         */
        double start = MPI_Wtime();
        PSpMat<OverlapPath>::MPI_DCCols N = Mult_AnXBn_DoubleBuff<MinPlusSR, OverlapPath, PSpMat<OverlapPath>::DCCols>(P, R);

        N.Prune(InvalidSRing(), true); /* GGGG: let's discuss this */

    #ifdef PDEBUG
        N.PrintInfo();  
    #endif

        P = N;
        timePR += MPI_Wtime() - start;

        OverlapPath id;

        /* GGGG: Ii is going to be true is the Ni dir entry corresponding to the Ri dir entry is non-zero, that is
            mark true edges that should be removed from R eventually because transitive.
            I would call this semiring something like "GreaterThenTR", but that's minor.
            */
        start = MPI_Wtime();
        PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, SpDCCols<int64_t, bool>>(F, N, TransitiveSelection(), false, id);

        I.Prune(BoolPrune(), true); /* GGGG: this is needed to remove entries in N that were smaller than F and thus resulted in an actual 0 in F */
        
    #ifdef PDEBUG
        I.PrintInfo();  
        tu.print_str("\n");
    #endif

        /* GGGG: make sure every transitive edge is correctly removed in the upper/lower triangular entry, too
            this would happen naturally on an error-free dataset 
            */

        /* I has some 1 entries, now IT has them, too */
        PSpMat<bool>::MPI_DCCols IT = I;

        IT.Transpose();

        if (!(IT == I)) /* symmetricize */
        {
            I += IT;
        }  

        timeI += MPI_Wtime() - start;
    #ifdef PDEBUG
        I.PrintInfo();  
    #endif

        start = MPI_Wtime();
        T += I; /* GGGG: add new transitive edges to T */

        cur = T.getnnz();
        timeT += MPI_Wtime() - start;

    #ifdef PDEBUG
        T.PrintInfo();  
        tu.print_str("\n");
    #endif

        /* GGGG: if nonzeros in T don't change for MAXITER iterations, then exit the loop 
            If in the test run this creates isses, e.g. it doesn't terminate, let's put a stricter exit condition */
    #ifdef BECONSERVATIVE
        if(cur == prev)
        {
            countidle++; 
        }
        else if(countidle > 0)
        {
            countidle = 0; // GGGG: reset the counter if cur != prev in this iteration but countidle > 0
        }
    #endif

        count++; // GGGG: this just keeps track of the total number of iterations but doesn't do anything about the termination condition

#ifdef BECONSERVATIVE
    } while (countidle < MAXITER);
#else
    } while (cur != prev);
#endif

    std::stringstream iss;
    iss << "TR took " << count << " iteration to complete!";
    tu.print_str(iss.str());

#ifdef PDEBUG
    T.PrintInfo();  
    T.ParallelWriteMM("transitive-matrix.mm", true);
#endif

    /* GGGG: this is not working as expected! there was a problem in the semiring but it's not fully fixed */
    R.PrintInfo();

    double start = MPI_Wtime();
    isLogicalNot = true; /* GGGG: I want the ones to be removed to be set to false, so I use the logical negation */ 
    bId = true; /* GGGG: this is critical, if this would be false, everything that would survive the EWIseApply would have dir == -1 as well and it's wrong
                    A non-transitive edge should keep its original direction! */
    R = EWiseApply<ReadOverlap, SpDCCols<int64_t, ReadOverlap>>(R, T, TransitiveRemoval(), isLogicalNot, bId);

#ifdef PDEBUG   
    tu.print_str("Matrix S, i.e. AAt post transitive reduction---BEFORE InvalidSRing Prune: ");
    R.PrintInfo();

    R.ParallelWriteMM("string-graph-before-prune.mm", true, ReadOverlapExtraHandler());
#endif

    R.Prune(InvalidSRing(), true);
    timeTR += MPI_Wtime() - start;

    tu.print_str("Matrix S, i.e. AAt post transitive reduction---AFTER InvalidSRing Prune: ");
    R.PrintInfo();

 #ifdef TIMINGTR
    R.ParallelWriteMM("string-matrix-after-tr.mm", true, ReadOverlapExtraHandler());
    double maxtimePR, maxtimeI, maxtimeT, maxtimeTR;

    MPI_Reduce(&timePR, &maxtimePR, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeI,  &maxtimeI,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeT,  &maxtimeT,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeTR, &maxtimeTR, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    std::stringstream tiss;
    tiss << "TransitiveReduction:TimePR (matmul) = "        << maxtimePR << "\n";
    tiss << "TransitiveReduction:TimeI  (element-wise) = "  <<  maxtimeI << "\n";
    tiss << "TransitiveReduction:TimeT  (element-wise) = "  <<  maxtimeT << "\n";
    tiss << "TransitiveReduction:TimeTR (element-wise) = "  << maxtimeTR << "\n";
    tu.print_str(tiss.str());
 #endif
 
}

#endif
