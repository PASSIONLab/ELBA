// Created by Saliya Ekanayake on 2019-02-17.

#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include <seqan/sequence.h>
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

  bool is_ready();

  void wait();

private:
  ushort k;
  ushort overlap;
  FastaData *fd = nullptr;

  uint64_t l_seq_count;
  uint64_t *l_seq_counts = nullptr;

  bool is_diagonal_cell = false;

  /*!
   * Global sequence count, which may be different from the total input sequence
   * count because some might get removed if their lengths are less than the
   * k-mer length.
   */
  uint64_t g_seq_count;
  uint64_t *g_seq_offsets = nullptr;

  std::vector<NbrData> my_nbrs;

  uint64_t row_seq_start_idx;
  uint64_t row_seq_end_idx;
  uint64_t col_seq_start_idx;
  uint64_t col_seq_end_idx;

  std::vector<seqan::Peptide *> row_seqs;
  std::vector<seqan::Peptide *> col_seqs;

  /*! recv counts and buffers */
  int recv_nbrs_count;
  std::vector<int> recv_nbrs_idxs;
  char **recv_nbrs_buffs = nullptr;
  uint64_t *recv_nbrs_buff_lengths = nullptr;
  MPI_Request *recv_nbrs_buffs_reqs = nullptr;
  MPI_Status *recv_nbrs_buffs_stats = nullptr;
  FastaData **recv_fds = nullptr;


  /*! send counts and buffers */
  int to_nbrs_count;
  MPI_Request *to_nbrs_buffs_reqs = nullptr;
  MPI_Status *to_nbrs_buffs_stat = nullptr;

  bool ready = false;

  std::shared_ptr<ParallelOps> parops;

  void collect_grid_seqs();

  void find_nbrs(int grid_rc_procs_count,
                 int grid_rc_id,
                 uint64_t avg_l_seq_count,
                 uint64_t avg_grid_seq_count,
                 uint64_t rc_seq_start_idx,
                 ushort rc_flag,
                 std::vector<NbrData> &my_nbrs);

  void push_seqs(int rc_flag, FastaData *fd, uint64_t seqs_count,
                 uint64_t seq_start_idx);

};


#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
