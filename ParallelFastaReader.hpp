// Created by Saliya Ekanayake on 1/6/19.

#ifndef LBL_DAL_PARALLEL_FASTA_READER_HPP
#define LBL_DAL_PARALLEL_FASTA_READER_HPP


#include "ParallelOps.hpp"
#include "DistributedFastaData.hpp"
#include "DebugUtils.hpp"

struct NbrData {
  ushort rc_flag;
  uint64_t nbr_rank;
  uint64_t nbr_seq_start_idx;
  uint64_t nbr_seq_end_idx;

  NbrData() = default;

  NbrData(ushort rc_flag, uint64_t nbr_rank, uint64_t nbr_seq_start_idx,
          uint64_t nbr_seq_end_idx) : rc_flag(rc_flag), nbr_rank(nbr_rank),
                                      nbr_seq_start_idx(nbr_seq_start_idx),
                                      nbr_seq_end_idx(nbr_seq_end_idx) {

  }
};

class ParallelFastaReader {
public:
  static std::unique_ptr<DistributedFastaData>
  read_and_distribute_fasta(const char *file, ushort overlap, ushort k,
                            const std::shared_ptr<ParallelOps> &parops);

private:
  static std::pair<char *, char *>
  find_grid_seqs(uint64_t g_seq_count, const uint64_t *g_seq_offsets,
                 uint64_t l_seq_count, const uint64_t *l_seq_counts,
                 const std::shared_ptr<ParallelOps> &parops);

  static void find_nbrs(
    uint64_t grid_rc_procs_count,
    uint64_t grid_rc_id,
    uint64_t g_seq_count,
    uint64_t avg_l_seq_count,
    uint64_t avg_grid_seq_count,
    uint64_t rc_seq_start_idx,
    ushort rc_flag,
    const uint64_t *g_seq_offsets,
    const uint64_t *l_seq_counts,
    const std::shared_ptr<ParallelOps> &parops,
    std::vector<NbrData> &nbrs
  );
};


#endif //LBL_DAL_PARALLEL_FASTA_READER_HPP
