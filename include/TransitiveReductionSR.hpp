// Created by Giulia Guidi on 09/04/20.

#ifndef __TER_DEFS_H__
#define __TER_DEFS_H__

#include "../include/Constants.hpp"
#include "../include/ParallelOps.hpp"
#include "../include/ParallelFastaReader.hpp"
#include "../include/Alphabet.hpp"
#include "../include/Utils.hpp"
#include "../include/DistributedPairwiseRunner.hpp"
#include "CombBLAS/CombBLAS.h"
#include "../include/cxxopts.hpp"
#include "../include/pw/SeedExtendXdrop.hpp"
#include "seqan/score/score_matrix_data.h"
#include "../include/pw/OverlapFinder.hpp"
#include "../include/pw/FullAligner.hpp"
#include "../include/pw/BandedAligner.hpp"
#include "../include/kmer/KmerOps.hpp"
#include "../include/kmer/KmerIntersectSR.hpp"

#include <sys/time.h> 
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>

#include <string>
#include <sstream>

using namespace std;
using namespace combblas;

#define FUZZ (10)
#define DEBUG

/*! GGGG: make these definitions consistent with the main code

typedef uint32_t dibella::CommonKmers;

// encoded int 
dibella::CommonKmers compose(const dibella::CommonKmers& suffix, const dibella::CommonKmers& dir) { return suffix << 2 | dir; }
// extract edge value 
dibella::CommonKmers val(const dibella::CommonKmers& me) { return me >> 2; }
// extract edge direction 
dibella::CommonKmers dir(const dibella::CommonKmers& me) { return me  & 3; }

dibella::CommonKmers min(const dibella::CommonKmers& arg1, const dibella::CommonKmers& arg2)
{
    if(val(arg2) < val(arg1)) return arg2;
    else return arg1;
}

dibella::CommonKmers max(const dibella::CommonKmers& arg1, const dibella::CommonKmers& arg2)
{
    if(val(arg2) > val(arg1)) return arg2;
    else return arg1;
}

template <class NT>
class PSpMat
{ 
public: 
	typedef SpDCCols <int64_t, NT> DCCols;
	typedef SpParMat <int64_t, NT, DCCols> MPI_DCCols;
};
*/

template <class T1, class T2, class OUT>
struct Bind2ndBiSRing : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const
    {
        return static_cast<OUT>(y);
    }
};

template <class T1, class T2, class OUT>
struct ReduceMBiSRing : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const
    {
        if(val(y) > val(x)) return static_cast<OUT>(y);
        else return static_cast<OUT>(x);
    }
};

template <class T, class OUT>
struct PlusFBiSRing : unary_function <T, OUT>
{
    OUT operator() (const T& x) const
    {
        return static_cast<OUT>(compose(val(x) + FUZZ, dir(x)));
    }
};

template <class T1, class T2, class OUT>
struct MinPlusBiSRing
{
	static OUT id() 			{ return std::numeric_limits<OUT>::max(); };
	static bool returnedSAID() 	{ return false; 	}
	static MPI_Op mpi_op() 		{ return MPI_MIN; 	};

	static OUT add(const OUT & arg1, const OUT & arg2)
	{
		return min(arg1, arg2);
	}
	static OUT multiply(const T1& arg1, const T2& arg2)
	{
        OUT res;
        if((dir(arg1) & 1) != (dir(arg2) & (1 << 1)))
        {
            OUT res = inf_plus<dibella::CommonKmers>(static_cast<dibella::CommonKmers>(val(arg1)), static_cast<dibella::CommonKmers>(val(arg2)));
            return compose(res, dir(arg2));
        } 
        else return id();
	}
	static void axpy(T1 a, const T2 & x, OUT & y)
	{
		y = min(y, multiply(a, x));
	}
};

template <class T1, class T2>
struct GreaterBinaryOp : binary_function <T1, T2, bool>
{
    bool operator() (const T1& x, const T2& y) const
    {
        if(val(x) >= val(y) && dir(y) == dir(x)) return true;
        else return false;
    }
};

template <class T1, class T2, class OUT>
struct MultiplyBinaryOp : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const { return static_cast<OUT>(compose(val(x) * val(y), dir(x))); }
};

template <class T>
struct ZeroUnaryOp : unary_function <T, bool>
{
    bool operator() (const T& x) const { if(x == 0) return true; else return false; }
};

#endif