// Created by Giulia Guidi on 04/02/21.

#ifndef __DIBELLA_SR_HPP__
#define __DIBELLA_SR_HPP__

#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"

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

// Bind length of 2nd but keep dir of 1st
template <class T1, class T2, class OUT>
struct Bind2ndBiSRing : binary_function <T1, T2, OUT>
{
    OUT operator() (const T1& x, const T2& y) const
    {
        OUT z;
        return static_cast<OUT>(compose(z, length(y) + FUZZ, dir(x)));
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

template <class T1, class T2, class OUT>
struct KeepShorterstSR : binary_function <T1, T2, OUT>
{
    // 
    OUT operator() (const T1& x, const T2& y) const
    {
        if(length(y) < length(x)) return static_cast<OUT>(y);
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

template <class T, class OUT>
struct IsEndOfContigSR : unary_function <T, OUT>
{
    OUT operator() (T& x) const
    {
        // return static_cast<OUT>(compose(x, length(x) + FUZZ, dir(x)));
    }
};

void tobinary(ushort n, int* arr) 
{ 
    int nbit = 2;
    for(int i = 0; i < nbit; i++)
    { 
        arr[i] = n % 2; 
        n = n / 2; 
    }
}

bool testdir(ushort dir1, ushort dir2, ushort& dir)
{
    ushort rbit, lbit;
    ushort start, end;

    int mybin1[2] = {0, 0};
    int mybin2[2] = {0, 0};

    if(dir1 != 0) tobinary(dir1, mybin1);
    if(dir2 != 0) tobinary(dir2, mybin2);

    rbit = mybin1[0];
    lbit = mybin2[1];

    if(rbit != lbit)
    {
        start = mybin1[1]; 
        end   = mybin2[0]; 

        if(start == 0)
        {
            if(end == 0) dir = 0;
            else dir = 1;
        }
        else
        {
            if(end == 0) dir = 2;
            else dir = 3;      
        }
        return true;
    }
    else return false;
}

dibella::CommonKmers concatenateseq(const dibella::CommonKmers& lhs, const dibella::CommonKmers& rhs, ushort& dir)
{
    dibella::CommonKmers res;

    // Giulia keep an eye on this; esay buggy; these sequences need to be on the same strand when concatenating them
    if(dir == 1 || dir == 0) 
    {
        char dest[strlen(lhs.seq)+strlen(rhs.seq)];

        strcpy(dest, lhs.seq);
        strcat(dest, rhs.seq);

        res.seq = dest;
    }
    else
    {
        char dest[strlen(lhs.seq)+strlen(rhs.seq)];

        strcpy(dest, rhs.seq);
        strcat(dest, lhs.seq);
        
        res.seq = dest;
    }

    return compose(res, length(lhs) + length(rhs), dir);
}

template <class T1, class T2, class OUT>
struct MinPlusBiSRing
{
	static OUT id() 			{ return std::numeric_limits<OUT>::max(); };
	static bool returnedSAID() 	{ return false; 	}
	static MPI_Op mpi_op() 		{ return MPI_MIN; 	};

	static OUT add(const OUT& arg1, const OUT& arg2)
	{
        if(arg1 < arg2) return arg1;
        else return arg2;
	}
	static OUT multiply(const T1& arg1, const T2& arg2)
	{
        OUT res;
        ushort mydir;

        if(testdir(dir(arg1), dir(arg2), mydir))
        {
            uint len = infplus(arg1, arg2);
            return compose(res, len, mydir);
        } 
        else return id();
	}
	static void axpy(T1 a, const T2& x, OUT& y)
	{   
		y = add(y, multiply(a, x));
	}
};

template <class T1, class T2, class OUT>
struct ContigSRing
{
	static OUT  id() 			{ 
        OUT res;
        return compose(res, 0, 0);
    };
	static bool returnedSAID() 	{ return false; 	}
	static MPI_Op mpi_op() 		{ return MPI_MIN; 	};

	static OUT add(const OUT & arg1, const OUT & arg2)
	{
        if(arg1 < arg2) return arg1;
        else return arg2;
	}
	static OUT multiply(const T1& arg1, const T2& arg2)
	{
        OUT res;
        ushort mydir;

        // printf("%d\t%d\n---\n", length(arg1), length(arg2));

        if(testdir(dir(arg1), dir(arg2), mydir))
        {
            return concatenateseq(arg1, arg2, mydir);
        } 
        else return id();
	}
	static void axpy(T1 a, const T2 & x, OUT & y)
	{   
		y = add(y, multiply(a, x));
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

template <class T1, class T2>
struct EqualBinaryOp : binary_function <T1, T2, bool>
{
    bool operator() (const T1& x, const T2& y) const
    {
        if(x == y) return true;
        else return false;
    }
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

template <class T, class OUT>
struct OverhangTSRing : unary_function <T, OUT>
{
    OUT operator() (const T& x) const
    {
        OUT xT = static_cast<OUT>(x);

        xT.overhang = x.overhangT;
        xT.overhangT = x.overhang;

        return xT;
    }
};

/*! Type definitions */
typedef MinPlusBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> MinPlusSR_t;
typedef ReduceMBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ReduceMSR_t;
typedef Bind2ndBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> Bind2ndSR_t;

#endif // __DIBELLA_SR_HPP__