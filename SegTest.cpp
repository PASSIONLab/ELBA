#include <iostream>
#include <random>
#include <chrono>
#include "CombBLAS/CombBLAS.h"

using namespace combblas;

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

struct CommonKmers {
  int count = 0;
  std::pair<int64_t, int64_t> first;
  std::pair<int64_t, int64_t> second;

  friend std::ostream &operator<<(std::ostream &os, const CommonKmers &m) {
    os << "|" << m.count << "(" << m.first.first << "," << m.first.second
       << ")(" <<
       m.second.first << "," << m.second.second << ")| ";
    return os;
  }
};


template<typename IN, typename OUT>
struct KmerIntersect {
  static OUT id() {
    OUT a;
    return a;
  }

  static bool returnedSAID() { return false; }

  static OUT add(const OUT &arg1, const OUT &arg2) {
    OUT res;
    res.count = arg1.count + arg2.count;
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

  static void axpy(IN a, const IN &x, OUT &y) {
    y = add(y, multiply(a, x));
  }

  static MPI_Op mpi_op() {
    static MPI_Op mpiop;
    static bool exists = false;
    if (exists)
      return mpiop;
    else {
      MPI_Op_create(MPI_func, true, &mpiop);
      exists = true;
      return mpiop;
    }
  }

  static void
  MPI_func(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
    for (int i = 0; i < *len; ++i) {
      *((OUT) inoutvec + i) = add(*((OUT) invec + i), *((OUT) inoutvec + 1));
    }

  }
};

int main(int argc, char **argv) {
  std::cout << "hello";
  int world_size, world_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//  assert(world_size == 4);

  /*!
   * For testing let's recreate matrices similar to the one that segfaults
   * Rank: 0 lrow_ids 146 lcol_ids 146 lvals 146
   * Rank: 1 lrow_ids 56 lcol_ids 56 lvals 56
   * Rank: 2 lrow_ids 100 lcol_ids 100 lvals 100
   * Rank: 3 lrow_ids 141 lcol_ids 141 lvals 141
   */

  unsigned int base = 20; // Represents the protein alphabet size
  unsigned int k = 2; // Kmer size
  unsigned long total_k = static_cast<unsigned int>(pow(base, k));

  // This is how many nnzs were there in each local matrix when segfaulted
  int nnzs[4] = {146, 56, 100, 141};

  std::vector<int64_t> lrow_ids(static_cast<unsigned long>(nnzs[world_rank]));
  std::vector<int64_t> lcol_ids(static_cast<unsigned long>(nnzs[world_rank]));
  std::vector<int64_t> lvals(static_cast<unsigned long>(nnzs[world_rank]));

  std::string fnames[4] = {"mat.0.txt", "mat.1.txt", "mat.2.txt", "mat.3.txt"};
  std::ifstream f(fnames[world_rank]);
  std::string v;
  int count = 0;
  while (f.good()){
    std::getline(f, v, ',');
    if (v.empty()) break;
    int val = std::stoi(v);
    lrow_ids[count] = val;
    std::getline(f, v, ',');
    val = std::stoi(v);
    lcol_ids[count] = val;
    std::getline(f, v);
    val = std::stoi(v);
    lvals[count] = val;
    ++count;
  }

  std::printf("Rank: %d lrow_ids %ld lcol_ids %ld lvals %ld\n",
              world_rank, lrow_ids.size(), lcol_ids.size(), lvals.size());

  std::shared_ptr<CommGrid> grid =
      std::make_shared<CommGrid>(MPI_COMM_WORLD, std::sqrt(world_size),std::sqrt(world_size));

  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, int64_t> dvals(lvals, grid);

  int m = 5, n = static_cast<int>(total_k);
  PSpMat<int64_t>::MPI_DCCols A(m, n, drows, dcols, dvals, false);

//  A.PrintInfo();

  auto At = A;
  At.Transpose();

  typedef KmerIntersect<int64_t, CommonKmers> KmerIntersectSR_t;

  PSpMat<CommonKmers>::MPI_DCCols C =
      Mult_AnXBn_Synch<KmerIntersectSR_t, CommonKmers, PSpMat<CommonKmers>::DCCols>(
          A, At);

  //seg faults here
  C.PrintInfo();

/* Hmm, someone else seems to call MPI_Finalize */
  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}

