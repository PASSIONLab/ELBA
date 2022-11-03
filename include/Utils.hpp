// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_SPARSEMAT_HPP
#define LBL_DAL_SPARSEMAT_HPP

#include <cstdint>
#include <utility>
#include <vector>
#include "CombBLAS/CombBLAS.h"
#include "Defines.hpp"
#include "kmer/Kmer.hpp"

// #define VERBOSE

using namespace std;

template <class T, std::size_t N>
ostream& operator<<(ostream& o, const array<T, N>& arr)
{
    copy(arr.cbegin(), arr.cend(), ostream_iterator<T>(o, " "));
    return o;
}

template <class NT>
class PSpMat
{
public:
  typedef combblas::SpTuples <int64_t, NT> Tuples;	
  typedef combblas::SpDCCols <int64_t, NT> DCCols;
  typedef combblas::SpParMat <int64_t, NT, DCCols> MPI_DCCols;
  typedef std::tuple<int64_t, int64_t, NT *> ref_tuples;
};

typedef vector<vector<array<char, 2>>>  VectorVectorChar;
typedef vector<vector<Kmer>>   VectorVectorKmer;
typedef vector<vector<ReadId>> VectorVectorReadId;
typedef vector<vector<PosInRead>> VectorVectorPos;
typedef vector<Kmer>  VectorKmer;
typedef vector<array<char, 2>> VectorChar;
typedef vector<PosInRead> VectorPos;
typedef vector<ReadId> VectorReadId;

#endif //LBL_DAL_SPARSEMAT_HPP
