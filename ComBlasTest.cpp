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

struct CommonKmers {
  int count = 0;
  std::pair<int64_t, int64_t> first;
  std::pair<int64_t, int64_t> second;

  friend std::ostream & operator << (std::ostream &os, const CommonKmers &m){
    os << "|"<<m.count<<"(" << m.first.first <<","<<m.first.second<< ")("<<
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
    res.first.first = arg1.first.first;
    res.first.second = arg1.first.second;
    res.second.first = arg2.first.first;
    res.second.second = arg2.first.second;
    return res;
  }
  static OUT multiply(const IN &arg1, const IN &arg2) {
    OUT a;
    a.count++;
    a.first.first = arg1;
    a.first.second = arg2;
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

int main(int argc, char** argv){
  std::cout<<"hello";
  int world_size, world_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // construct objects
  MPI_Comm world = MPI_COMM_WORLD;


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
  /*switch (world_rank){
    case 0:{
      lrow_ids.push_back(0); lcol_ids.push_back(0); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(5); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(6); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(1); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(2); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(3); lvals.push_back(1);
      lrow_ids.push_back(1); lcol_ids.push_back(6); lvals.push_back(1);
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
  }*/


  /*
   * For testing let's create a sparse matrix as follows
   *
   * p0 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 0| 1 x x 4 x 6 7
   * 1| x 9 10 11 x x 14
   * 2| x x x x 19 x 21
   * -------------
   *
   * p1 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 3| x x x 25 26 27 28
   * 4| x 30 31 32 33 x x
   * 5| 36 37 38 39 40 41 x
   *
   * p2 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 6| x x x x x 48 49
   * 7| x x 52 x 54 x 56
   *
   * p3 finds k-mers for this part
   *    0 1 2 3 4 5 6
   * -------------
   * 8| 57 58 59 x x 62 x
   * 9| 64 x x x x 68 69
   *
   */
  switch (world_rank){
    case 0:{
      lrow_ids.push_back(0); lcol_ids.push_back(0); lvals.push_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(3); lvals.push_back(4);
      lrow_ids.push_back(0); lcol_ids.push_back(5); lvals.push_back(6);
      lrow_ids.push_back(0); lcol_ids.push_back(6); lvals.push_back(7);
      lrow_ids.push_back(1); lcol_ids.push_back(1); lvals.push_back(9);
      lrow_ids.push_back(1); lcol_ids.push_back(2); lvals.push_back(10);
      lrow_ids.push_back(1); lcol_ids.push_back(3); lvals.push_back(11);
      lrow_ids.push_back(1); lcol_ids.push_back(6); lvals.push_back(14);
      lrow_ids.push_back(2); lcol_ids.push_back(4); lvals.push_back(19);
      lrow_ids.push_back(2); lcol_ids.push_back(6); lvals.push_back(21);
      break;
    }
    case 1:{
      lrow_ids.push_back(3); lcol_ids.push_back(3); lvals.push_back(25);
      lrow_ids.push_back(3); lcol_ids.push_back(4); lvals.push_back(26);
      lrow_ids.push_back(3); lcol_ids.push_back(5); lvals.push_back(27);
      lrow_ids.push_back(3); lcol_ids.push_back(6); lvals.push_back(28);
      lrow_ids.push_back(4); lcol_ids.push_back(1); lvals.push_back(30);
      lrow_ids.push_back(4); lcol_ids.push_back(2); lvals.push_back(31);
      lrow_ids.push_back(4); lcol_ids.push_back(3); lvals.push_back(32);
      lrow_ids.push_back(4); lcol_ids.push_back(4); lvals.push_back(33);
      lrow_ids.push_back(5); lcol_ids.push_back(0); lvals.push_back(36);
      lrow_ids.push_back(5); lcol_ids.push_back(1); lvals.push_back(37);
      lrow_ids.push_back(5); lcol_ids.push_back(2); lvals.push_back(38);
      lrow_ids.push_back(5); lcol_ids.push_back(3); lvals.push_back(39);
      lrow_ids.push_back(5); lcol_ids.push_back(4); lvals.push_back(40);
      lrow_ids.push_back(5); lcol_ids.push_back(5); lvals.push_back(41);
      break;
    }
    case 2:{
      lrow_ids.push_back(6); lcol_ids.push_back(5); lvals.push_back(48);
      lrow_ids.push_back(6); lcol_ids.push_back(6); lvals.push_back(49);
      lrow_ids.push_back(7); lcol_ids.push_back(2); lvals.push_back(52);
      lrow_ids.push_back(7); lcol_ids.push_back(4); lvals.push_back(54);
      lrow_ids.push_back(7); lcol_ids.push_back(6); lvals.push_back(56);
      break;
    }
    case 3:{
      lrow_ids.push_back(8); lcol_ids.push_back(0); lvals.push_back(57);
      lrow_ids.push_back(8); lcol_ids.push_back(1); lvals.push_back(58);
      lrow_ids.push_back(8); lcol_ids.push_back(2); lvals.push_back(59);
      lrow_ids.push_back(8); lcol_ids.push_back(5); lvals.push_back(62);
      lrow_ids.push_back(9); lcol_ids.push_back(0); lvals.push_back(64);
      lrow_ids.push_back(9); lcol_ids.push_back(5); lvals.push_back(69);
      lrow_ids.push_back(9); lcol_ids.push_back(6); lvals.push_back(70);
      break;
    }
  }

  std::shared_ptr<CommGrid> grid = std::make_shared<CommGrid>(MPI_COMM_WORLD, 2, 2);
  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, int64_t> dvals(lvals, grid);

  int m = 10, n = 7;
  PSpMat<int64_t >::MPI_DCCols A(m, n, drows, dcols, dvals, false);

//  A.PrintInfo();

  /*std::printf("World rank |%d| \n", world_rank);
  int64_t nnz = A.getnnz();
  if (world_rank == 0){
    std::printf("\nwow nnz %lli\n", nnz);
  }*/

  auto At = A;
  At.Transpose();

//  PSpMat<int64_t>::MPI_DCCols C =
//      Mult_AnXBn_Synch<PlusTimesSRing<int64_t, int64_t >, int64_t, PSpMat<int64_t >::DCCols>(A, At);


  typedef KmerIntersect<int64_t, CommonKmers> KmerIntersectSR_t;

  PSpMat<CommonKmers>::MPI_DCCols C =
      Mult_AnXBn_Synch<KmerIntersectSR_t, CommonKmers, PSpMat<CommonKmers>::DCCols>(A, At);
  C.PrintInfo();

  //
  /*PSpMat<CommonKmers>::DCCols* arrays = C.seqptr();
  auto a = arrays->GetArrays();
  auto idx_arrays = a.indarrs;
  std::cout<<idx_arrays.size();*/







/* Hmm, someone else seems to call MPI_Finalize */
//  MPI_Finalize();
}
