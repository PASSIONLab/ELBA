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

  /*
   * For testing let's create a sparse matrix as follows
   *
   * p0 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 0| 1 x x 1 x 1 1
   * 1| x 1 1 1 x x 1
   * 2| x x x x 1 x 1
   * -------------
   *
   * p1 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 3| x x x 1 1 1 1
   * 4| x 1 1 1 1 x x
   * 5| 1 1 1 1 1 1 x
   *
   * p2 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 6| x x x x x 1 1
   * 7| x x 1 x 1 x 1
   *
   * p3 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 8| 1 1 1 x x 1 x
   * 9| 1 x x x x 1 1
   *
   */

  std::vector<int64_t> lrow_ids, lcol_ids, lvals;

  assert(world_size == 4);
  switch (world_rank){
    case 0:{
      lrow_ids.push_back(0); lcol_ids.push_back(0); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(5); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(6); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(1); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(2); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(2); lcol_ids.push_back(6); lvals.push_back(1);
      lrow_ids.push_back(2); lcol_ids.push_back(4); lvals.push_back(1);
      lrow_ids.push_back(2); lcol_ids.push_back(6); lvals.push_back(1);
      break;
    }
    case 1:{
      lrow_ids.push_back(3); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(3); lcol_ids.push_back(4); lvals.push_back(1);
      lrow_ids.push_back(3); lcol_ids.push_back(5); lvals.push_back(1);
      lrow_ids.push_back(3); lcol_ids.push_back(6); lvals.push_back(1);
      lrow_ids.push_back(4); lcol_ids.push_back(1); lvals.push_back(1);
      lrow_ids.push_back(4); lcol_ids.push_back(2); lvals.push_back(1);
      lrow_ids.push_back(4); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(4); lcol_ids.push_back(4); lvals.push_back(1);
      lrow_ids.push_back(5); lcol_ids.push_back(0); lvals.push_back(1);
      lrow_ids.push_back(5); lcol_ids.push_back(1); lvals.push_back(1);
      lrow_ids.push_back(5); lcol_ids.push_back(2); lvals.push_back(1);
      lrow_ids.push_back(5); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(5); lcol_ids.push_back(4); lvals.push_back(1);
      lrow_ids.push_back(5); lcol_ids.push_back(5); lvals.push_back(1);
      break;
    }
    case 2:{
      lrow_ids.push_back(6); lcol_ids.push_back(5); lvals.push_back(1);
      lrow_ids.push_back(6); lcol_ids.push_back(6); lvals.push_back(1);
      lrow_ids.push_back(7); lcol_ids.push_back(2); lvals.push_back(1);
      lrow_ids.push_back(7); lcol_ids.push_back(4); lvals.push_back(1);
      lrow_ids.push_back(7); lcol_ids.push_back(6); lvals.push_back(1);
      break;
    }
    case 3:{
      lrow_ids.push_back(8); lcol_ids.push_back(0); lvals.push_back(1);
      lrow_ids.push_back(8); lcol_ids.push_back(1); lvals.push_back(1);
      lrow_ids.push_back(8); lcol_ids.push_back(2); lvals.push_back(1);
      lrow_ids.push_back(8); lcol_ids.push_back(5); lvals.push_back(1);
      lrow_ids.push_back(9); lcol_ids.push_back(0); lvals.push_back(1);
      lrow_ids.push_back(9); lcol_ids.push_back(5); lvals.push_back(1);
      lrow_ids.push_back(9); lcol_ids.push_back(6); lvals.push_back(1);
      break;
    }
  }

  std::shared_ptr<CommGrid> grid = std::make_shared<CommGrid>(MPI_COMM_WORLD, 2, 2);
  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);


/* Hmm, someone else seems to call MPI_Finalize */
//  MPI_Finalize();
}
