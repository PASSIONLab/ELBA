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

std::unique_ptr<CT<Overlap>::PSpParMat> TransitiveReduction(CT<Overlap>::PSpParMat R);

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

Overlap opmin(const Overlap& e1, const Overlap& e2);

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

#endif
