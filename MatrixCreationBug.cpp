// Created by Saliya Ekanayake on 10/16/19.
#include <iostream>
#include "CombBLAS/CombBLAS.h"

using namespace combblas;

struct ME{
  short cost;
  unsigned short offset;

  ME():cost(0), offset(0){}
  ME(short cost, unsigned short offset):cost(cost), offset(offset){}

  friend std::ostream& operator<<(std::ostream& os, const ME& me){
    os << "(" << me.cost << ", " << me.offset << ")";
    return os;
  }

  ME operator+(const ME& me) const{
    ME me2(cost,offset);
    return me2;
  }

  bool operator<(const ME& me) const{
    return true;
  }
};


/* Copied from MultTest.cpp in CombBLAS */
// Simple helper class for declarations: Just the numerical type is templated
// The index type and the sequential matrix type stays the same for the whole code
// In this case, they are "int" and "SpDCCols"
template<class NT>
class PSpMat {
public:
  typedef SpDCCols<int64_t, NT> DCCols;
  typedef SpParMat<int64_t, NT, DCCols> MPI_DCCols;
};

int64_t pack(short a, unsigned short b)
{
  return (int64_t)((((uint64_t)a)<<16)+(uint32_t)b);
}

short unpackA(int x)
{
  return (short)(((uint64_t)x)>>16);
}

unsigned short unpackB(int x)
{
  return (unsigned short)(((uint64_t)x)&0xffff);
}

void test_gen_A_Pack(int world_rank, int world_size){
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
  std::vector<int64_t> lrow_ids, lcol_ids;
  std::vector<int64_t> lvals;
  assert(world_size == 4);

  switch (world_rank){
    case 0:{
      lrow_ids.push_back(0); lcol_ids.push_back(0); lvals.emplace_back(pack(1,1));
      lrow_ids.push_back(0); lcol_ids.push_back(3); lvals.emplace_back(pack(4, 4));
      lrow_ids.push_back(0); lcol_ids.push_back(5); lvals.emplace_back(pack(6, 6));
      lrow_ids.push_back(0); lcol_ids.push_back(6); lvals.emplace_back(pack(7, 7));
      lrow_ids.push_back(1); lcol_ids.push_back(1); lvals.emplace_back(pack(9, 9));
      lrow_ids.push_back(1); lcol_ids.push_back(2); lvals.emplace_back(pack(10, 10));
      lrow_ids.push_back(1); lcol_ids.push_back(3); lvals.emplace_back(pack(11, 11));
      lrow_ids.push_back(1); lcol_ids.push_back(6); lvals.emplace_back(pack(14, 14));
      lrow_ids.push_back(2); lcol_ids.push_back(4); lvals.emplace_back(pack(1, 19));
      lrow_ids.push_back(2); lcol_ids.push_back(6); lvals.emplace_back(pack(1, 21));
      break;
    }
    case 1:{
      lrow_ids.push_back(3); lcol_ids.push_back(3); lvals.emplace_back(pack(5, 25));
      lrow_ids.push_back(3); lcol_ids.push_back(4); lvals.emplace_back(pack(6, 26));
      lrow_ids.push_back(3); lcol_ids.push_back(5); lvals.emplace_back(pack(7, 27));
      lrow_ids.push_back(3); lcol_ids.push_back(6); lvals.emplace_back(pack(8, 28));
      lrow_ids.push_back(4); lcol_ids.push_back(1); lvals.emplace_back(pack(0, 30));
      lrow_ids.push_back(4); lcol_ids.push_back(2); lvals.emplace_back(pack(1, 31));
      lrow_ids.push_back(4); lcol_ids.push_back(3); lvals.emplace_back(pack(2, 32));
      lrow_ids.push_back(4); lcol_ids.push_back(4); lvals.emplace_back(pack(3, 33));
      lrow_ids.push_back(5); lcol_ids.push_back(0); lvals.emplace_back(pack(6, 36));
      lrow_ids.push_back(5); lcol_ids.push_back(1); lvals.emplace_back(pack(7, 37));
      lrow_ids.push_back(5); lcol_ids.push_back(2); lvals.emplace_back(pack(8, 38));
      lrow_ids.push_back(5); lcol_ids.push_back(3); lvals.emplace_back(pack(9, 39));
      lrow_ids.push_back(5); lcol_ids.push_back(4); lvals.emplace_back(pack(0, 40));
      lrow_ids.push_back(5); lcol_ids.push_back(5); lvals.emplace_back(pack(1, 41));
      break;
    }
    case 2:{
      lrow_ids.push_back(6); lcol_ids.push_back(5); lvals.emplace_back(pack(8, 48));
      lrow_ids.push_back(6); lcol_ids.push_back(6); lvals.emplace_back(pack(9, 49));
      lrow_ids.push_back(7); lcol_ids.push_back(2); lvals.emplace_back(pack(2, 52));
      lrow_ids.push_back(7); lcol_ids.push_back(4); lvals.emplace_back(pack(4, 54));
      lrow_ids.push_back(7); lcol_ids.push_back(6); lvals.emplace_back(pack(6, 56));
      break;
    }
    case 3:{
      lrow_ids.push_back(8); lcol_ids.push_back(0); lvals.emplace_back(pack(7, 57));
      lrow_ids.push_back(8); lcol_ids.push_back(1); lvals.emplace_back(pack(8, 58));
      lrow_ids.push_back(8); lcol_ids.push_back(2); lvals.emplace_back(pack(9, 59));
      lrow_ids.push_back(8); lcol_ids.push_back(5); lvals.emplace_back(pack(2, 62));
      lrow_ids.push_back(9); lcol_ids.push_back(0); lvals.emplace_back(pack(4, 64));
      lrow_ids.push_back(9); lcol_ids.push_back(5); lvals.emplace_back(pack(9, 69));
      lrow_ids.push_back(9); lcol_ids.push_back(6); lvals.emplace_back(pack(0, 70));
      break;
    }
  }

  std::shared_ptr<CommGrid> grid
      = std::make_shared<CommGrid>(MPI_COMM_WORLD, 2, 2);

  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, int64_t> dvals(lvals, grid);

  int m = 10, n = 7;
  PSpMat<int64_t>::MPI_DCCols A(m, n, drows, dcols, dvals, false);

  A.PrintInfo();
}

void test_gen_A_short(int world_rank, int world_size){
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
  std::vector<int64_t> lrow_ids, lcol_ids;
  std::vector<short> lvals;
  assert(world_size == 4);

  switch (world_rank){
    case 0:{
      lrow_ids.push_back(0); lcol_ids.push_back(0); lvals.emplace_back(1);
      lrow_ids.push_back(0); lcol_ids.push_back(3); lvals.emplace_back(4);
      lrow_ids.push_back(0); lcol_ids.push_back(5); lvals.emplace_back(6);
      lrow_ids.push_back(0); lcol_ids.push_back(6); lvals.emplace_back(7);
      lrow_ids.push_back(1); lcol_ids.push_back(1); lvals.emplace_back(9);
      lrow_ids.push_back(1); lcol_ids.push_back(2); lvals.emplace_back(10);
      lrow_ids.push_back(1); lcol_ids.push_back(3); lvals.emplace_back(11);
      lrow_ids.push_back(1); lcol_ids.push_back(6); lvals.emplace_back(14);
      lrow_ids.push_back(2); lcol_ids.push_back(4); lvals.emplace_back(1);
      lrow_ids.push_back(2); lcol_ids.push_back(6); lvals.emplace_back(1);
      break;
    }
    case 1:{
      lrow_ids.push_back(3); lcol_ids.push_back(3); lvals.emplace_back(5);
      lrow_ids.push_back(3); lcol_ids.push_back(4); lvals.emplace_back(6);
      lrow_ids.push_back(3); lcol_ids.push_back(5); lvals.emplace_back(7);
      lrow_ids.push_back(3); lcol_ids.push_back(6); lvals.emplace_back(8);
      lrow_ids.push_back(4); lcol_ids.push_back(1); lvals.emplace_back(0);
      lrow_ids.push_back(4); lcol_ids.push_back(2); lvals.emplace_back(1);
      lrow_ids.push_back(4); lcol_ids.push_back(3); lvals.emplace_back(2);
      lrow_ids.push_back(4); lcol_ids.push_back(4); lvals.emplace_back(3);
      lrow_ids.push_back(5); lcol_ids.push_back(0); lvals.emplace_back(6);
      lrow_ids.push_back(5); lcol_ids.push_back(1); lvals.emplace_back(7);
      lrow_ids.push_back(5); lcol_ids.push_back(2); lvals.emplace_back(8);
      lrow_ids.push_back(5); lcol_ids.push_back(3); lvals.emplace_back(9);
      lrow_ids.push_back(5); lcol_ids.push_back(4); lvals.emplace_back(0);
      lrow_ids.push_back(5); lcol_ids.push_back(5); lvals.emplace_back(1);
      break;
    }
    case 2:{
      lrow_ids.push_back(6); lcol_ids.push_back(5); lvals.emplace_back(8);
      lrow_ids.push_back(6); lcol_ids.push_back(6); lvals.emplace_back(9);
      lrow_ids.push_back(7); lcol_ids.push_back(2); lvals.emplace_back(2);
      lrow_ids.push_back(7); lcol_ids.push_back(4); lvals.emplace_back(4);
      lrow_ids.push_back(7); lcol_ids.push_back(6); lvals.emplace_back(6);
      break;
    }
    case 3:{
      lrow_ids.push_back(8); lcol_ids.push_back(0); lvals.emplace_back(7);
      lrow_ids.push_back(8); lcol_ids.push_back(1); lvals.emplace_back(8);
      lrow_ids.push_back(8); lcol_ids.push_back(2); lvals.emplace_back(9);
      lrow_ids.push_back(8); lcol_ids.push_back(5); lvals.emplace_back(2);
      lrow_ids.push_back(9); lcol_ids.push_back(0); lvals.emplace_back(4);
      lrow_ids.push_back(9); lcol_ids.push_back(5); lvals.emplace_back(9);
      lrow_ids.push_back(9); lcol_ids.push_back(6); lvals.emplace_back(0);
      break;
    }
  }

  std::shared_ptr<CommGrid> grid
      = std::make_shared<CommGrid>(MPI_COMM_WORLD, 2, 2);

  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, short> dvals(lvals, grid);

  int m = 10, n = 7;
  PSpMat<short>::MPI_DCCols A(m, n, drows, dcols, dvals, false);

  A.PrintInfo();
}


void test_gen_A_ME(int world_rank, int world_size){
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
  std::vector<int64_t> lrow_ids, lcol_ids;
  std::vector<ME> lvals;
  assert(world_size == 4);

  switch (world_rank){
    case 0:{
      lrow_ids.push_back(0); lcol_ids.push_back(0); lvals.emplace_back(1, 1);
      lrow_ids.push_back(0); lcol_ids.push_back(3); lvals.emplace_back(4, 4);
      lrow_ids.push_back(0); lcol_ids.push_back(5); lvals.emplace_back(6, 6);
      lrow_ids.push_back(0); lcol_ids.push_back(6); lvals.emplace_back(7, 7);
      lrow_ids.push_back(1); lcol_ids.push_back(1); lvals.emplace_back(9, 9);
      lrow_ids.push_back(1); lcol_ids.push_back(2); lvals.emplace_back(10, 10);
      lrow_ids.push_back(1); lcol_ids.push_back(3); lvals.emplace_back(11, 11);
      lrow_ids.push_back(1); lcol_ids.push_back(6); lvals.emplace_back(14, 14);
      lrow_ids.push_back(2); lcol_ids.push_back(4); lvals.emplace_back(1, 19);
      lrow_ids.push_back(2); lcol_ids.push_back(6); lvals.emplace_back(1, 21);
      break;
    }
    case 1:{
      lrow_ids.push_back(3); lcol_ids.push_back(3); lvals.emplace_back(5, 25);
      lrow_ids.push_back(3); lcol_ids.push_back(4); lvals.emplace_back(6, 26);
      lrow_ids.push_back(3); lcol_ids.push_back(5); lvals.emplace_back(7, 27);
      lrow_ids.push_back(3); lcol_ids.push_back(6); lvals.emplace_back(8, 28);
      lrow_ids.push_back(4); lcol_ids.push_back(1); lvals.emplace_back(0, 30);
      lrow_ids.push_back(4); lcol_ids.push_back(2); lvals.emplace_back(1, 31);
      lrow_ids.push_back(4); lcol_ids.push_back(3); lvals.emplace_back(2, 32);
      lrow_ids.push_back(4); lcol_ids.push_back(4); lvals.emplace_back(3, 33);
      lrow_ids.push_back(5); lcol_ids.push_back(0); lvals.emplace_back(6, 36);
      lrow_ids.push_back(5); lcol_ids.push_back(1); lvals.emplace_back(7, 37);
      lrow_ids.push_back(5); lcol_ids.push_back(2); lvals.emplace_back(8, 38);
      lrow_ids.push_back(5); lcol_ids.push_back(3); lvals.emplace_back(9, 39);
      lrow_ids.push_back(5); lcol_ids.push_back(4); lvals.emplace_back(0, 40);
      lrow_ids.push_back(5); lcol_ids.push_back(5); lvals.emplace_back(1, 41);
      break;
    }
    case 2:{
      lrow_ids.push_back(6); lcol_ids.push_back(5); lvals.emplace_back(8, 48);
      lrow_ids.push_back(6); lcol_ids.push_back(6); lvals.emplace_back(9, 49);
      lrow_ids.push_back(7); lcol_ids.push_back(2); lvals.emplace_back(2, 52);
      lrow_ids.push_back(7); lcol_ids.push_back(4); lvals.emplace_back(4, 54);
      lrow_ids.push_back(7); lcol_ids.push_back(6); lvals.emplace_back(6, 56);
      break;
    }
    case 3:{
      lrow_ids.push_back(8); lcol_ids.push_back(0); lvals.emplace_back(7, 57);
      lrow_ids.push_back(8); lcol_ids.push_back(1); lvals.emplace_back(8, 58);
      lrow_ids.push_back(8); lcol_ids.push_back(2); lvals.emplace_back(9, 59);
      lrow_ids.push_back(8); lcol_ids.push_back(5); lvals.emplace_back(2, 62);
      lrow_ids.push_back(9); lcol_ids.push_back(0); lvals.emplace_back(4, 64);
      lrow_ids.push_back(9); lcol_ids.push_back(5); lvals.emplace_back(9, 69);
      lrow_ids.push_back(9); lcol_ids.push_back(6); lvals.emplace_back(0, 70);
      break;
    }
  }

  std::shared_ptr<CommGrid> grid
  = std::make_shared<CommGrid>(MPI_COMM_WORLD, 2, 2);

  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, ME> dvals(lvals, grid);

  int m = 10, n = 7;
  PSpMat<ME>::MPI_DCCols A(m, n, drows, dcols, dvals, false);

  A.PrintInfo();
}

void test_gen_A(int world_rank, int world_size){
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
  std::vector<int64_t> lrow_ids, lcol_ids, lvals;
  assert(world_size == 4);

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

  std::shared_ptr<CommGrid> grid = std::make_shared<CommGrid>(MPI_COMM_WORLD, 2,
                                                              2);
  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, int64_t> dvals(lvals, grid);

  int m = 10, n = 7;
  PSpMat<int64_t>::MPI_DCCols A(m, n, drows, dcols, dvals, false);

  A.PrintInfo();
}

int main(int argc, char **argv) {
  int world_size, world_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  // construct objects
  MPI_Comm world = MPI_COMM_WORLD;



//  test_gen_A(world_rank, world_size);
  test_gen_A_ME(world_rank, world_size);
//  test_gen_A_Pack(world_rank, world_size);
//  test_gen_A_short(world_rank, world_size);
  MPI_Finalize();
}
