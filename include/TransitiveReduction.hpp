
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

#ifndef FUZZ
#define FUZZ 1000
#endif

int intplus(int a, int b)
{
    return (a == MAX_INT || b == MAX_INT) ? MAX_INT : a + b;
}

struct InvalidSRing : unary_function<ReadOverlap, ReadOverlap>
{
    bool operator() (const ReadOverlap& x) { return (x.valid == 0 || x.direction() == -1); }
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

        xT.refix(1);

        return xT;
    }
};

ReadOverlap omin(const ReadOverlap& e1, const ReadOverlap& e2)
{
    ReadOverlap e = ReadOverlap();

    for (int i = 0; i < 4; ++i)
        e.sfx[i] = std::min(e1.sfx[i], e2.sfx[i]);
    
    return e;
}

struct MinPlusSR
{
    static ReadOverlap id() { return ReadOverlap(); }
    static bool returnedSAID() { return false; } /* what does this do? */
    static MPI_Op mpi_op() { return MPI_MIN; }   /* what does this do? */

    static ReadOverlap add(const ReadOverlap& e1, const ReadOverlap& e2)
    {
        return omin(e1, e2);
    }

    static ReadOverlap multiply(const ReadOverlap& e1, const ReadOverlap& e2)
    {
        ReadOverlap e = ReadOverlap();

        int indove = e1.direction()&1;

        e.sfx[0] = intplus(e1.sfx[indove],   e2.sfx[2*(!indove)]  ); /* >--@ @--< */
        e.sfx[1] = intplus(e1.sfx[indove],   e2.sfx[2*(!indove)+1]); /* >--@ @--> */
        e.sfx[2] = intplus(e1.sfx[indove+2], e2.sfx[2*(!indove)]  ); /* <--@ @--< */
        e.sfx[3] = intplus(e1.sfx[indove+2], e2.sfx[2*(!indove)+1]); /* <--@ @--> */

        return e;
    }

    static void axpy(ReadOverlap a, const ReadOverlap& x, ReadOverlap& y)
    {
        y = omin(y, multiply(a, x));
    }
};

struct PlusFuzzSRing : unary_function<ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& x) const
    {
        ReadOverlap fuzzed = x;

        for (int i = 0; i < 4; ++i)
            fuzzed.sfx[i] = intplus(fuzzed.sfx[i], FUZZ);

        return fuzzed;
    }
};

struct TransitiveSelection : binary_function<ReadOverlap, ReadOverlap, bool>
{
    bool operator() (const ReadOverlap& x, const ReadOverlap& y) const
    {
        int xdir = x.direction();
        return (x.sfx[xdir] >= y.sfx[xdir]);
    }
};

struct TransitiveRemoval : binary_function<ReadOverlap, bool, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& x, const bool& y)
    {
        if (!y) x.valid = 0;
        return x;
    }
};

void TransitiveReduction(SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>& R, TraceUtils tu)
{
    SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> RT = R;
    RT.Transpose();
    RT.Apply(TransposeSRing());

    if (!(RT == R)) R += RT;

    R.Prune(InvalidSRing(), true);

    R.ParallelWriteMM("R.mm", true, ReadOverlapMMHandler());

    int nnz, prev;
    do {

        tu.print_str("Matrix R, i.e. AAt, mid transitive reduction: ");
        R.PrintInfo();

        prev = R.getnnz();

        SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> Rc = R;
        SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> N = Mult_AnXBn_DoubleBuff<MinPlusSR, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>(R, Rc);

        N.Prune(InvalidSRing(), true);

        SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> M = R;

        M.Apply(PlusFuzzSRing());

        ReadOverlap id;

        SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> I = EWiseApply<bool, SpDCCols<int64_t, bool>>(M, N, TransitiveSelection(), false, id);

        R = EWiseApply<ReadOverlap, SpDCCols<int64_t, ReadOverlap>>(R, I, TransitiveRemoval(), true, true);

        R.Prune(InvalidSRing(), true);

        nnz = R.getnnz();
    } while (nnz != prev);

    R.ParallelWriteMM("S.mm", true, ReadOverlapMMHandler());

    tu.print_str("Matrix R, i.e. AAt after transitive reduction: ");
    R.PrintInfo();
}

#endif
