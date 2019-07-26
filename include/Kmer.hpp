// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_KMER_HPP
#define LBL_DAL_KMER_HPP

#include <cstdint>
#include <vector>
#include <cmath>
#include "ParallelOps.hpp"
#include "Alphabet.hpp"
#include "Types.hpp"
#include "ScoreMat.hpp"


struct CommonKmers {
  /*! The number of common kmers between two sequences.
   * The maximum could be floor((l-k)/s)+1, where
   * l is the sequence length, k is the kmer length, and
   * s is the stride. Since l is within 2^16-1 (unsigned short max)
   * we can represent the count as unsigned short as well.
   */
  ushort count;

  /*! The position within the sequence, which is
   * much less than 2^16 - 1 for proteins
   */
  std::pair<ushort, ushort> first;
  std::pair<ushort, ushort> second;

  CommonKmers() : count(1) {
  }

  explicit CommonKmers(ushort count) : count(count){
  }

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
    OUT res(arg1.count + arg2.count);
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
                        uint64_t end_offset_inclusive, ushort k, ushort s,
                        bool add_substitute_kmers,
                        std::map<uint64_t, std::vector<uint64_t>*> &kmer_to_subs_kmers,
                        pisa::ScoreMatrix &score_mat, Alphabet &alp,
                        uvec_64 &lcol_ids, uvec_16 &lvals,
                        const std::shared_ptr<ParallelOps> parops) {

  auto num_kmers = static_cast<ushort>((floor((len - k) * 1.0 / s)) + 1);
  // TODO: Saliya - this can be improved using bit operators
  ushort base = alp.size;
  uint64_t kcode = 0;
  int count = 0;
  char cap_c;
  for (uint64_t i = start_offset; i <= ((end_offset_inclusive - k) + 1); i += s) {
    kcode = 0;
    if (!add_substitute_kmers && count == num_kmers) break;
    for (uint64_t j = i; j < i + k; ++j) {
      /*! Efficient than using pow() */
      cap_c = *(seq + j);
      if (cap_c > 96 && cap_c < 123){
        // small case character, so make it uppercase.
        cap_c = cap_c - 32;
      }
      kcode = kcode * base + alp.char_to_code[cap_c];
    }

    if (!add_substitute_kmers) {
      ++count;
      lcol_ids.push_back(kcode);
      /*! Offset is relative to the sequence start, so unsigned short is
       * good enough.
       */
      lvals.push_back(static_cast<ushort &&>(i - start_offset));
    } else {
      if (kmer_to_subs_kmers.find(kcode) == kmer_to_subs_kmers.end()) {
        /*! Encountering this kmer for the first time, so let's find its
         * substitute kmers and add them.
         */

        uint64_t start_idx = 0;
        auto *subs_kcodes = new std::vector<uint64_t>();
        subs_kcodes->push_back(0);
        for (uint64_t j = i; j < i + k; ++j) {
          /*! Efficient than using pow() */
          cap_c = *(seq + j);
          if (cap_c > 96 && cap_c < 123) {
            // small case character, so make it uppercase.
            cap_c = cap_c - 32;
          }

          uint64_t subs_kcodes_size = subs_kcodes->size();
          for (int m = start_idx; m < subs_kcodes_size; ++m){
            uint64_t subs_kcode = (*subs_kcodes)[m];
            for (auto const &sub : (*score_mat.base_to_subtitutes[cap_c])){
              subs_kcodes->push_back(subs_kcode*base+alp.char_to_code[sub]);
            }
          }
          start_idx = subs_kcodes_size;
        }

        subs_kcodes->erase(subs_kcodes->begin(), subs_kcodes->begin()+start_idx);
        kmer_to_subs_kmers[kcode] = subs_kcodes;
      }

      /*! This kmer was previously encountered or just discovered now in the
       * last step. Either way we already know what
       * are its substitute kmers. Let's add them.
       */
      std::vector<uint64_t> *subs_kmers = kmer_to_subs_kmers[kcode];
      for (auto const  &subs_kcode : (*subs_kmers)){
        lcol_ids.push_back(subs_kcode);
        lvals.push_back(static_cast<ushort &&>(i - start_offset));
        ++count;
      }

    }
  }

  if (!add_substitute_kmers && count != num_kmers) {
    fprintf(stderr,
            "ERROR: rank: %d, count:%d numk: %d len: %d k: %d s: %d soff: %llu eoffinc: %llu\n",
            parops->world_proc_rank, count, num_kmers, len, k, s,
            start_offset, end_offset_inclusive);
    fflush(stderr);
  }

  return count;
}


#endif //LBL_DAL_KMER_HPP
