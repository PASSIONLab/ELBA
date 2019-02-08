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
  typedef combblas::SpDCCols <uint64_t, NT> DCCols;
  typedef combblas::SpParMat <uint64_t, NT, DCCols> MPI_DCCols;
};

#endif //LBL_DAL_SPARSEMAT_HPP
