// Created by Giulia Guidi on 09/04/20.

#ifndef __TER_DEFS_H__
#define __TER_DEFS_H__

#include "../include/kmer/CommonKmers.hpp"
#include "../include/Utils.hpp"

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

/*! Namespace declarations */
using namespace combblas;
using namespace std;

#define FUZZ (1000)
#define DEBUG

/** Given a biridrected graph, an edge v ?-? x can only be considered transitive given a pair of edges v ?-? w ?-? x if:
 * (1) The two heads adjacent to w have opposite orientation: 
 *      2nd bit != 1st bit in MinPlus semiring B = A^2 such as 01 and 01 or 10 and 10;
 * (2) The heads adjacent to v in v ?-? w and v ?-? x have the same orientation, and
 * (3) The heads adjacent to x in v ?-? x and w ?-? x have the same orientation:
 *      1st and 2nd bit in M == 1st and 2nd bit in B during I = M >= B.
*/

dibella::CommonKmers compose(dibella::CommonKmers& me, const uint& suffix, const ushort& dir)
{ 
    me.overhang = suffix << 2 | dir;
    return me; 
}

void tobinary(ushort n, int* arr) 
{ 
    int nbit = 2;
    for(int i = 0; i < nbit; i++)
    { 
        arr[i] = n % 2; 
        n = n / 2; 
    }
} 

ushort updatedir(ushort dir1, ushort dir2)
{
    ushort res;

    if((dir1 == 0 & (dir2 == 1 || dir2 == 2)) 
        || ((dir1 == 1 || dir1 == 2) & dir2 == 0))
    {
        res = 0;
    }
    else if((dir1 == 3 & (dir2 == 1 || dir2 == 2)) 
            || ((dir1 == 1 || dir1 == 2) & dir2 == 3))
    {
        res = 3;      
    }
    else if((dir1 == 0 & dir2 == 3) 
            || (dir1 == 1 & (dir2 == 1 || dir2 == 2)))
    {
        res = 1;     
    }
    else if((dir1 == 3 & dir2 == 0) 
            || (dir1 == 2 & (dir2 == 1 || dir2 == 2)))
    {
        res = 2;     
    }

    return res;
}

uint length(const dibella::CommonKmers& me) { return me.overhang >> 2; }
ushort  dir(const dibella::CommonKmers& me) { return me.overhang  & 3; }

dibella::CommonKmers min(const dibella::CommonKmers& arg1, const dibella::CommonKmers& arg2) {
    if(length(arg2) < length(arg1)) return arg2;
    else return arg1;
}

dibella::CommonKmers max(const dibella::CommonKmers& arg1, const dibella::CommonKmers& arg2) {
    if(length(arg2) > length(arg1)) return arg2;
    else return arg1;
}

const uint infplus(const dibella::CommonKmers& a, const dibella::CommonKmers& b) {
	uint inf = std::numeric_limits<uint>::max();
    if (length(a) == inf || length(b) == inf) {
    	return inf;
    }
    return length(a) + length(b);
}

// GGGG: best (minimum overhang) overlap semiring on vector
template <class T1, class T2, class OUT>
struct BestOverlapVSRing : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const
    {
        if(length(y) < length(x)) return static_cast<OUT>(y);
        else return static_cast<OUT>(x);
    }
};

// GGGG: best (minimum overhang) overlap semiring on matrix
template <class T1, class T2, class OUT>
struct BestOverlapMSRing : binary_function <T1, T2, OUT>
{
    OUT operator() (T1& x, const T2& y) const
    {
        /* I want to only keep entry corresponding to the best overlap for that column */
        if(length(x) != length(y)) x.overhang = 0;
        return static_cast<OUT>(x);
    }
};

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
        if(length(y) > length(x)) return static_cast<OUT>(y);
        else return static_cast<OUT>(x);
    }
};

template <class T, class OUT>
struct PlusFBiSRing : unary_function <T, OUT>
{
    OUT operator() (T& x) const
    {
        return static_cast<OUT>(compose(x, length(x) + FUZZ, dir(x)));
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

        // if 01 and 10 are flexible then the only case we don't compute is if they are both 00 or both 11
        // the rest need to be computed and direction updated
        // this is gonna generate extra entries but technically the invalid one shouldn't affect the output
        // if we find out that they do, we should implement global consistency using broadcasting (see slides)
        if((dir(arg1) != 0 & dir(arg2) != 0) 
            || (dir(arg1) != 3 & dir(arg2) != 3))
        {
            uint   len   = infplus(arg1, arg2);
            ushort mydir = updatedir(dir(arg1), dir(arg2));

            return compose(res, len, mydir);
        } 
        else 
        {
            return id();
        }
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
        if(length(x) >= length(y) && dir(y) == dir(x)) return true;
        else return false;
    }
};

template <class T1, class T2, class OUT>
struct MultiplyBinaryOp : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const { return static_cast<OUT>(compose(length(x) * length(y), dir(x))); }
};

template <class T>
struct ZeroUnaryOp : unary_function <T, bool>
{
    bool operator() (const T& x) const { if(x == 0) return true; else return false; }
};

template <class T, class T2>
struct EWiseMulOp : binary_function <T, T2, T>
{
    T operator() (T& x, const T2& y) const
    {
        if(!y) x.overhang = 0;
        return x;
    }
};

template <class T>
struct ZeroOverhangSR : unary_function <T, bool>
{
    bool operator() (const T& x) const { if(x.overhang == 0) return true; else return false; }
};

/*! Type definitions */
typedef MinPlusBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> MinPlusSR_t;
typedef ReduceMBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ReduceMSR_t;
typedef Bind2ndBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> Bind2ndSR_t;

/* TR main function */
void TransitiveReduction(PSpMat<dibella::CommonKmers>::MPI_DCCols& B, TraceUtils tu)
{
    PSpMat<dibella::CommonKmers>::MPI_DCCols BT = B;
    BT.Transpose();

    if(!(BT == B))
    {
        // @@GGGG-TODO: write overload function for dibella::CommonKmers
        B += BT;
        // B.Apply();
    }
    B.PrintInfo();

    uint nnz, prev;
    double timeA2 = 0, timeC = 0, timeI = 0, timeA = 0;

    /* Gonna iterate on B until there are no more transitive edges to remove */
    do
    {
      prev = B.getnnz();

      /* Find two-hops neighbors
       * C = B^2
       */
      double start = MPI_Wtime();
      PSpMat<dibella::CommonKmers>::MPI_DCCols F = B;
      PSpMat<dibella::CommonKmers>::MPI_DCCols C = Mult_AnXBn_DoubleBuff<MinPlusSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, F);
      timeA2 += MPI_Wtime() - start;
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix C = B^2: ");
      C.PrintInfo();
  #endif
    
      start = MPI_Wtime();
      FullyDistVec<int64_t, dibella::CommonKmers> vA(B.getcommgrid());

      dibella::CommonKmers id; 
      vA = B.Reduce(Row, ReduceMSR_t(), id);
      vA.Apply(PlusFBiSRing<dibella::CommonKmers, dibella::CommonKmers>());

      F.DimApply(Row, vA, Bind2ndSR_t());
      timeC += MPI_Wtime() - start;
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix F = B + FUZZ: ");
      F.PrintInfo();
  #endif

      /* Find transitive edges that can be removed
       * I = F >= C 
       */
      start = MPI_Wtime();
      bool isLogicalNot = false;
      PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(F, C, GreaterBinaryOp<dibella::CommonKmers, dibella::CommonKmers>(), isLogicalNot, id);

      /* Prune potential zero-valued nonzeros */
      I.Prune(ZeroUnaryOp<bool>(), true);
      timeI += MPI_Wtime() - start;
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix I = F >= B: ");
      I.PrintInfo();
  #endif

      /* Remove transitive edges
       * B = B .* not(I)
       */ 
      start = MPI_Wtime();
      isLogicalNot = true;
      B = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, I, EWiseMulOp<dibella::CommonKmers, bool>(), isLogicalNot, true);

      /* Prune zero-valued overhang */
      B.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);
      timeA += MPI_Wtime() - start;

  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix B = B .* not(I): ");
      B.PrintInfo();
  #endif
      nnz = B.getnnz();
       
    } while (nnz != prev);

    tu.print_str("Matrix B, i.e AAt after transitive reduction: ");
    B.PrintInfo();

 #ifdef DIBELLA_DEBUG
    double maxtimeA2, maxtimeC, maxtimeI, maxtimeA;
    
    MPI_Reduce(&timeA2, &maxtimeA2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeC,  &maxtimeC,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeI,  &maxtimeI,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeA,  &maxtimeA,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(myrank == 0)
    {
      std::cout << "TransitiveReduction:TimeA2 = " << maxtimeA2 << std::endl;
      std::cout << "TransitiveReduction:TimeC  = " <<  maxtimeC << std::endl;
      std::cout << "TransitiveReduction:TimeI  = " <<  maxtimeI << std::endl;
      std::cout << "TransitiveReduction:TimeA  = " <<  maxtimeA << std::endl;
    }
 #endif
}

#endif