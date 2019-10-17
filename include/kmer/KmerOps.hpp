// Created by Saliya Ekanayake on 10/15/19.

#ifndef LBL_PISA_KMEROPS_HPP
#define LBL_PISA_KMEROPS_HPP

#include <unordered_set>
#include "MatrixEntry.hpp"
#include "../Types.hpp"
#include "../ParallelOps.hpp"
#include "../Alphabet.hpp"
#include "../ScoreMat.hpp"
#include "../Utils.hpp"
#include "../DistributedFastaData.hpp"

namespace pisa {
  class KmerOps {
  public:
    static PSpMat<MatrixEntry>::MPI_DCCols generate_A(uint64_t seq_count,
        std::shared_ptr<DistributedFastaData>& dfd, ushort k, ushort s,
        Alphabet &alph, const std::shared_ptr<ParallelOps>& parops,
        const std::shared_ptr<TimePod>& tp);

    static uint64_t add_kmers(const char *seq, ushort len, uint64_t start_offset,
                              uint64_t end_offset_inclusive, ushort k, ushort s,
                              Alphabet &alp, uvec_64 &lcol_ids, std::vector<MatrixEntry> &lvals,
                              const std::shared_ptr<ParallelOps>& parops) {

      auto num_kmers = static_cast<ushort>((floor((len - k) * 1.0 / s)) + 1);
      // TODO: Saliya - this can be improved using bit operators
      ushort base = alp.size;
      uint64_t kcode = 0;
      uint64_t count = 0;
      char cap_c;
      std::unordered_set<uint64_t> kmers_in_sequence;
      for (uint64_t i = start_offset;
           i <= ((end_offset_inclusive - k) + 1); i += s) {
        kcode = 0;
        if (count == num_kmers) break;
        for (uint64_t j = i; j < i + k; ++j) {
          /*! Efficient than using pow() */
          cap_c = *(seq + j);
          if (cap_c > 96 && cap_c < 123) {
            // small case character, so make it uppercase.
            cap_c = cap_c - 32;
          }
          kcode = kcode * base + alp.char_to_code[cap_c];
        }

        ++count;
        lcol_ids.push_back(kcode);
        /*! Offset is relative to the sequence start, so unsigned short is
         * good enough. */
        lvals.emplace_back(0, static_cast<ushort &&>(i - start_offset));
      }

      if (count != num_kmers) {
        fprintf(stderr,
                "ERROR: kmerop: add_kmers(): rank: %d, count:%d numk: %d len: %d k: %d s: %d soff: %llu eoffinc: %llu\n",
                parops->world_proc_rank, count, num_kmers, len, k, s,
                start_offset, end_offset_inclusive);
        fflush(stderr);
      }
      return count;
    }
  };
}


#endif //LBL_PISA_KMEROPS_HPP
