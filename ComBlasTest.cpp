//
// Created by Saliya Ekanayake on 12/21/18.
//

#include <iostream>
#include "CombBLAS/CombBLAS.h"

using namespace combblas;

/* Copied from MultTest.cpp in CombBLAS */
// Simple helper class for declarations: Just the numerical type is templated
// The index type and the sequential matrix type stays the same for the whole code
// In this case, they are "int" and "SpDCCols"
template <class NT>
class PSpMat
{
public:
  typedef SpDCCols < int64_t, NT > DCCols;
  typedef SpParMat < int64_t, NT, DCCols > MPI_DCCols;
};

int main(int argc, char** argv){
  std::cout<<"hello";
  int world_size, world_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // construct objects
  MPI_Comm world = MPI_COMM_WORLD;
  PSpMat<double>::MPI_DCCols A(world);

/* Hmm, someone else seems to call MPI_Finalize */
//  MPI_Finalize();
}
