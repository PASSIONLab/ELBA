// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_KMER_HPP
#define LBL_DAL_KMER_HPP

#include <cstdint>
#include <vector>
#include <cmath>
#include "Alphabet.hpp"

typedef std::vector<int64_t> vec64_t;

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

int add_kmers(const char *seq, int len, int start_offset,
              int enf_offset_inclusive, int k, int s, Alphabet &alp,
              vec64_t &lcol_ids, vec64_t &lvals){
  auto num_kmers = static_cast<int>(floor((len - k)*1.0 / s)) + 1;
  // TODO: Saliya - this can be improved using bit operators
  int64_t base = alp.size;
  int64_t kcode = 0;
  int count = 0;
  for (int i = start_offset; i <= enf_offset_inclusive; i+=s){
    kcode = 0;
    if (count == num_kmers) break;
    for (int j = i; j < i+k; ++j){
      /*! Efficient than using pow() */
      kcode = kcode * base + alp.char_to_code[*(seq + j)];
    }
    ++count;
    lcol_ids.push_back(kcode);
    lvals.push_back(i);
  }

  return num_kmers;
}



#endif //LBL_DAL_KMER_HPP
