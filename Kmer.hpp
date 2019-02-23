// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_KMER_HPP
#define LBL_DAL_KMER_HPP

#include <cstdint>
#include <vector>
#include <cmath>
#include "ParallelOps.hpp"
#include "Alphabet.hpp"
#include "Types.hpp"


struct CommonKmers {
  /*! The number of common kmers between two sequences.
   * The maximum could be floor((l-k)/s)+1, where
   * l is the sequence length, k is the kmer length, and
   * s is the stride. Since l is within 2^16-1 (unsigned short max)
   * we can represent the count as unsigned short as well.
   */
  ushort count = 0;

  /*! The position within the sequence, which is
   * much less than 2^16 - 1 for proteins
   */
  std::pair<ushort, ushort> first;
  std::pair<ushort, ushort> second;

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

static ushort add_kmers(const char *seq, ushort len, uint64_t start_offset,
                        uint64_t enf_offset_inclusive, ushort k, ushort s,
                        Alphabet &alp,
                        uvec_64 &lcol_ids, uvec_16 &lvals,
                        const std::shared_ptr<ParallelOps> parops) {
  auto num_kmers = static_cast<ushort>((floor((len - k) * 1.0 / s)) + 1);
  // TODO: Saliya - this can be improved using bit operators
  ushort base = alp.size;
  uint64_t kcode = 0;
  int count = 0;
  for (uint64_t i = start_offset; i <= enf_offset_inclusive; i += s) {
    kcode = 0;
    if (count == num_kmers) break;
    for (uint64_t j = i; j < i + k; ++j) {
      /*! Efficient than using pow() */
      kcode = kcode * base + alp.char_to_code[*(seq + j)];
    }
    ++count;
    lcol_ids.push_back(kcode);
    /*! Offset is relative to the sequence start, so unsigned short is
     * good enough.
     */
    lvals.push_back(static_cast<ushort &&>(i - start_offset));
  }

  if (count != num_kmers) {
    fprintf(stderr,
            "ERROR: rank: %d, count:%d numk: %d len: %d k: %d s: %d soff: %llu eoffinc: %llu\n",
            parops->world_proc_rank, count, num_kmers, len, k, s,
            start_offset, enf_offset_inclusive);
    fflush(stderr);
  }

  return num_kmers;
}


#endif //LBL_DAL_KMER_HPP
