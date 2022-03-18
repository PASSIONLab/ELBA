
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
    bool operator() (const ReadOverlap& x) { return !x.isvalid(); }
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

struct PlusFuzzSRing : unary_function<ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& x) const
    {
        ReadOverlap fuzzed = x;
        fuzzed.sfx += FUZZ;

        return fuzzed;
    }
};

//struct FlipReverseCoordinates : unary_function<ReadOverlap, ReadOverlap>
//{
//    ReadOverlap operator() (ReadOverlap& x) const
//    {
//        if (!x.rc) return x;
//
//        int swap = x.b[1];
//        x.b[1] = x.l[1] - x.e[1];
//        x.e[1] = x.l[1] - swap;
//
//        return x;
//    }
//};

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

//struct OverlapAdd : unary_function<ReadOverlap, ReadOverlap>
//{
//    ReadOverlap operator() (ReadOverlap& x) const
//    {
//        ReadOverlap out = x;
//        out.sfx += static_cast<int64_t>(out.overlap * 0.05);
//        return out;
//    }
//};

void TransitiveReduction(SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>& R, TraceUtils tu)
{

    R.ParallelWriteMM("overlap-graph-tri.mm", true, ReadOverlapExtraHandler());

    SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> RT = R;
    RT.Transpose();
    RT.Apply(TransposeSRing());

    if (!(RT == R)) R += RT;

    R.Prune(InvalidSRing(), true);

    R.ParallelWriteMM("overlap-graph.mm", true, ReadOverlapExtraHandler());

    SpParMat<int64_t, OverlapPath, SpDCCols<int64_t, OverlapPath>> Nc = R;
    SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> T = R;
    T.Prune(ZeroPrune());

    int64_t prev, cur;

    std::stringstream ss;
    int i = 1;

    SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> M = R;
    M.Apply(PlusFuzzSRing());

    do {
        tu.print_str("Matrix R, i.e. AAt, mid transitive reduction: ");
        R.PrintInfo();

        ss.str("");
        ss << i;

        prev = T.getnnz();
        SpParMat<int64_t, OverlapPath, SpDCCols<int64_t, OverlapPath>> N = Mult_AnXBn_DoubleBuff<MinPlusSR, OverlapPath, SpDCCols<int64_t, OverlapPath>>(Nc, R);
        N.Prune(InvalidSRing(), true);
        Nc = N;

        OverlapPath id;
        SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> I = EWiseApply<bool, SpDCCols<int64_t, bool>>(M, N, TransitiveSelection(), false, id);
        SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> It = I;
        It.Transpose();

        I.EWiseMult(It, false);
        SpParMat<int64_t, bool, SpDCCols<int64_t, bool>> Tc = T;

        T = EWiseApply<bool, SpDCCols<int64_t, bool>>(Tc, I, [](bool x, bool y) { return !(x && y); }, true, false);
        cur = T.getnnz();

    } while (i++ <= 5);
    //} while (0);
    //} while (prev != cur);

    R = EWiseApply<ReadOverlap, SpDCCols<int64_t, ReadOverlap>>(R, T, TransitiveRemoval(), false, false);
    R.Prune(InvalidSRing());

    //R.ParallelWriteMM("S.mm", true, ReadOverlapHandler());
    R.ParallelWriteMM("string-graph.mm", true, ReadOverlapExtraHandler());

    tu.print_str("Matrix R, i.e. AAt after transitive reduction: ");
    R.PrintInfo();
}

#endif
