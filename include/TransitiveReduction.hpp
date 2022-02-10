
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
    if (a == MAX_INT || b == MAX_INT) return MAX_INT;
    return a + b;
}

ReadOverlap omin(const ReadOverlap& u, const ReadOverlap& v)
{
        ReadOverlap uxv = ReadOverlap();

        uxv.sfx[0] = std::min(u.sfx[0], v.sfx[0]);
        uxv.sfx[1] = std::min(u.sfx[1], v.sfx[1]);
        uxv.sfx[2] = std::min(u.sfx[2], v.sfx[2]);
        uxv.sfx[3] = std::min(u.sfx[3], v.sfx[3]);

        return uxv;
}

struct MinPlusSR
{
    static ReadOverlap id() { return ReadOverlap(); }
    static bool returnedSAID() { return false; } /* what does this do? */
    static MPI_Op mpi_op() { return MPI_MIN; }   /* what does this do? */

    static ReadOverlap add(const ReadOverlap& u, const ReadOverlap& v)
    {
        return omin(u, v);
    }

    static ReadOverlap multiply(const ReadOverlap& u, const ReadOverlap& v)
    {
        ReadOverlap uxv = ReadOverlap();

        int udir = u.direction();

        if (udir&1) {
            uxv.sfx[0] = intplus(u.sfx[1], v.sfx[0]); /* >------>  >------< */
            uxv.sfx[1] = intplus(u.sfx[1], v.sfx[1]); /* >------>  >------> */
            uxv.sfx[2] = intplus(u.sfx[3], v.sfx[0]); /* <------>  >------< */
            uxv.sfx[3] = intplus(u.sfx[3], v.sfx[1]); /* <------>  >------> */
        } else {
            uxv.sfx[0] = intplus(u.sfx[0], v.sfx[2]);
            uxv.sfx[1] = intplus(u.sfx[0], v.sfx[3]);
            uxv.sfx[2] = intplus(u.sfx[2], v.sfx[2]);
            uxv.sfx[3] = intplus(u.sfx[2], v.sfx[3]);
        }

        return uxv;
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

void TransitiveReduction(SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>& R)
{
    SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> RT = R;
    RT.Transpose();
    RT.Apply(TransposeSRing());

    if (!(RT == R)) R += RT;

    R.Prune(InvalidSRing(), true);

    R.ParallelWriteMM("R.mm", true, ReadOverlapHandler());

    int nnz, prev;
    do {

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
}


#endif
