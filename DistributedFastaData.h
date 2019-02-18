// Created by Saliya Ekanayake on 2019-02-17.

#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H


#include "FastaData.hpp"
#include "ParallelOps.hpp"

struct NbrData {
  ushort rc_flag;
  int owner_rank;
  int nbr_rank;
  uint64_t nbr_seq_start_idx;
  uint64_t nbr_seq_end_idx;

  NbrData() = default;

  NbrData(ushort rc_flag, int owner_rank, int nbr_rank,
          uint64_t nbr_seq_start_idx,
          uint64_t nbr_seq_end_idx)
    : rc_flag(rc_flag),
      owner_rank(owner_rank),
      nbr_rank(nbr_rank),
      nbr_seq_start_idx(nbr_seq_start_idx),
      nbr_seq_end_idx(nbr_seq_end_idx) {

  }
};

class DistributedFastaData {
public:
  ~DistributedFastaData();

  DistributedFastaData(const char *file, ushort overlap, ushort k,
                       const std::shared_ptr<ParallelOps> &parops);

  uint64_t global_count();
  bool grid_ready();
  void wait();

private:
  ushort k;
  ushort overlap;
  FastaData *fd = nullptr;

  uint64_t l_seq_count;
  uint64_t *l_seq_counts = nullptr;

  /*!
   * Global sequence count, which may be different from the total input sequence
   * count because some might get removed if their lengths are less than the
   * k-mer length.
   */
  uint64_t g_seq_count;
  uint64_t *g_seq_offsets = nullptr;

  int recv_nbrs_count;
  int to_nbrs_count;
//  uint64_t *recv_nbrs_buff_lengths = nullptr;
  char **recv_nbrs_buffs = nullptr;
  MPI_Request *recv_nbrs_buffs_reqs = nullptr;
  MPI_Status *recv_nbrs_buffs_stats = nullptr;
  MPI_Request *to_nbrs_buffs_reqs = nullptr;
  MPI_Status *to_nbrs_bufss_stat = nullptr;

  bool ready = false;

  std::shared_ptr<ParallelOps> parops;

  void collect_grid_seqs();
  void find_nbrs(const int grid_rc_procs_count,
                 const int grid_rc_id,
                 const uint64_t avg_l_seq_count,
                 const uint64_t avg_grid_seq_count,
                 const uint64_t rc_seq_start_idx,
                 const ushort rc_flag,
                 std::vector<NbrData> &my_nbrs);

};


#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
