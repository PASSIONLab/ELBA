#ifndef KMER_INTERSECT_H_
#define KMER_INTERSECT_H_

#include "ReadOverlap.hpp"
#include <mpi.h>
#include <cstdlib>
#include <cmath>
#include <functional>

struct KmerIntersect
{
    static ReadOverlap id()
    {
        ReadOverlap a;
        return a;
    }

    static bool returnedSAID() { return false; }

    static ReadOverlap add(const ReadOverlap& arg1, const ReadOverlap& arg2)
    {
        ReadOverlap res(arg1.count + arg2.count);

        res.begQs[0] = arg1.begQs[0];
        res.begQs[1] = arg2.begQs[0];

        res.begTs[0] = arg1.begTs[0];
        res.begTs[1] = arg2.begTs[0];

        return res;
    }

    static ReadOverlap multiply(const PosInRead& arg1, const PosInRead& arg2)
    {
        ReadOverlap a;

        a.begQs[0] = arg1;
        a.begTs[0] = arg2;

        return a;
    }

    static void axpy(PosInRead a, const PosInRead& x, ReadOverlap& y)
    {
        y = add(y, multiply(a, x));
    }

    // static MPI_Op mpi_op()
    // {
        // static MPI_Op mpiop;
        // static bool exists = false;

        // if (!exists)
        // {
            // MPI_Op_create(MPI_func, true, &mpiop);
            // exists = true;
        // }

        // return mpiop;
    // }

    // static void MPI_func(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype)
    // {
        // for (int i = 0; i < *len; ++i)
        // {
            // *((ReadOverlap*) inoutvec + i) = add(*((ReadOverlap*) invec + i), *((ReadOverlap*) inoutvec + 1));
        // }
    // }
};

#endif
