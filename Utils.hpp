// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_SPARSEMAT_HPP
#define LBL_DAL_SPARSEMAT_HPP

#include <cstdint>
#include <utility>
#include "CombBLAS/CombBLAS.h"

template <class NT>
class PSpMat
{
public:
  typedef combblas::SpDCCols < int64_t, NT > DCCols;
  typedef combblas::SpParMat < int64_t, NT, DCCols > MPI_DCCols;
};

struct CommonKmers {
  int count = 0;
  std::pair<int64_t, int64_t> first;
  std::pair<int64_t, int64_t> second;

  friend std::ostream & operator << (std::ostream &os, const CommonKmers &m){
    os << "|(" << m.first.first <<","<<m.first.second<< ")("<<
       m.second.first<<","<<m.second.second<<")| ";
    return os;
  }
};

template <typename IN, typename OUT>
struct KmerIntersect
{
  static OUT id(){
    OUT a;
    return a;
  }
  static bool returnedSAID() { return false; }

  static OUT add(const OUT& arg1, const OUT& arg2) {
    OUT res;
    res.count = arg1.count + arg2.count;
    // TODO: perhaps improve this late with something that'll check how far
    // apart are the kmers.
    res.second.first = arg2.first.first;
    res.second.second = arg2.first.second;
    return arg1;
  }
  static OUT multiply(const IN &arg1, const IN &arg2) {
    OUT a;
    a.count++;
    a.first.first = arg1;
    a.first.second = arg1;
    return a;
  }
  static void axpy(IN a, const IN & x, OUT & y)
  {
    y = add(y, multiply(a, x));
  }

  static MPI_Op mpi_op()
  {
    static MPI_Op mpiop;
    static bool exists = false;
    if (exists)
      return mpiop;
    else
    {
      MPI_Op_create(MPI_func, true, &mpiop);
      exists = true;
      return mpiop;
    }
  }

  static void MPI_func(void * invec, void * inoutvec, int * len, MPI_Datatype *datatype)
  {
    for (int i = 0; i < *len; ++i){
      *((OUT)inoutvec+i) = add(*((OUT)invec+i), *((OUT)inoutvec+1));
    }

  }
};

#endif //LBL_DAL_SPARSEMAT_HPP
