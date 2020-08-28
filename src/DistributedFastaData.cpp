// Created by Saliya Ekanayake on 2019-02-17.

#include "../include/DistributedFastaData.hpp"
#include "../include/ParallelFastaReader.hpp"


DistributedFastaData::~DistributedFastaData() {
  for (auto &row_seq : row_seqs) {
    delete (row_seq);
  }

  if (!is_diagonal_cell) {
    /*! If this was a diagonal cell then both the row_seqs and col_seqs
     * would point to the same sequences, so deleting row_seqs is enough.
     */
    for (auto &col_seq : col_seqs) {
      delete (col_seq);
    }
  }

  if (recv_fds != nullptr) {
    for (int i = 0; i < recv_nbrs_count; ++i) {
      delete (recv_fds[i]);
    }
    delete[]recv_fds;
  }
  delete[](to_nbrs_buffs_stat);
  delete[](to_nbrs_buffs_reqs);
  delete[](recv_nbrs_buffs_stats);
  delete[](recv_nbrs_buffs_reqs);
  if (recv_nbrs_buffs != nullptr) {
    for (int i = 0; i < recv_nbrs_count; ++i) {
      if (recv_nbrs_buffs[i] != nullptr) {
        delete[](recv_nbrs_buffs[i]);
      }
    }
    delete[](recv_nbrs_buffs);
  }
  delete[](recv_nbrs_buff_lengths);
  delete (fd->buffer());
  delete (fd);
  delete[](l_seq_counts);
  delete[](g_seq_offsets);
}

DistributedFastaData::DistributedFastaData(
  const char *file, const char* idx_map_file,
  uint64_t overlap, ushort k,
  const std::shared_ptr<ParallelOps> &parops,
  const std::shared_ptr<TimePod> &tp, TraceUtils tu)
  : overlap(overlap), k(k), parops(parops), tp(tp), tu(tu) {

  is_diagonal_cell =
    parops->grid->GetRankInProcRow() == parops->grid->GetRankInProcCol();

  tp->times["StartDfd:PfrReadFasta()"] = std::chrono::system_clock::now();
  char *buff;
  uint64_t l_start, l_end;
  ParallelFastaReader::read_fasta(file, overlap, parops->world_proc_rank,
                                  parops->world_procs_count, buff, l_start,
                                  l_end);

  tp->times["EndDfd:PfrReadFasta()"] = std::chrono::system_clock::now();

  tp->times["StartDfd:newFD()"] = tp->times["EndDfd:PfrReadFasta()"];
  fd = new FastaData(buff, k, l_start, l_end, tp, tu);
  l_seq_count = fd->local_count();

  tp->times["EndDfd:newFD()"] = std::chrono::system_clock::now();

#ifndef NDEBUG
  {
    std::string title = "l_seq_count";
    std::string msg = std::to_string(l_seq_count);
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

  l_seq_counts = new uint64_t[parops->world_procs_count];
  MPI_Allgather(&l_seq_count, 1, MPI_UINT64_T, l_seq_counts,
                1, MPI_UINT64_T, MPI_COMM_WORLD);

  g_seq_offsets = new uint64_t[parops->world_procs_count];
  g_seq_offsets[0] = 0;
  for (int i = 1; i < parops->world_procs_count; ++i) {
    g_seq_offsets[i] = g_seq_offsets[i - 1] + l_seq_counts[i - 1];
  }

  /*! Total sequence count. This could be different from the original
   * input sequence count as we remove sequences less than k-mer length.
   */
  g_seq_count = g_seq_offsets[parops->world_procs_count - 1] +
                l_seq_counts[parops->world_procs_count - 1];

  uint64_t orig_l_seq_count = fd->orig_local_count();
  orig_g_seq_offset = 0;
  MPI_Exscan(&orig_l_seq_count, &orig_g_seq_offset, 1,
          MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);

  write_idx_map(idx_map_file);

#ifndef NDEBUG
  std::printf("Rank: %d orig_l_seq_count: %lld\n",
          parops->world_proc_rank, orig_l_seq_count);
#endif
  collect_grid_seqs();
}

void DistributedFastaData::write_idx_map(const char* idx_map_file){
    uint64_t g_seq_offset = g_seq_offsets[parops->world_proc_rank];
    uvec_64 *del_idxs = fd->deleted_indices();
    uint64_t del_count = del_idxs->size();
    uint64_t cur_dels = 0;
    uint64_t orig_idx;
    std::stringstream ss;
    for (uint64_t i = 0; i < l_seq_count; ++i){
        if (cur_dels < del_count && i+cur_dels == (*del_idxs)[cur_dels]){
            // This index has been deleted originally.
            ++cur_dels;
        }

        orig_idx = i + cur_dels + orig_g_seq_offset;
        ss << (i + g_seq_offset) << "," << orig_idx <<"\n";
    }
    parops->write_file_in_parallel(idx_map_file, ss.str());
}

void DistributedFastaData::collect_grid_seqs() {
  int pr = parops->grid->GetGridRows();
  int pc = parops->grid->GetGridCols();

  /*! It might seem wrong to use rank_in_col for row_id and rank_in_row for
   * col_id but it is correct. These are ranks when the column and row worlds
   * are viewed separately.
   */
  int grid_row_id = parops->grid->GetRankInProcCol();
  int grid_col_id = parops->grid->GetRankInProcRow();

  uint64_t avg_l_seq_count = g_seq_count / parops->world_procs_count;
  uint64_t avg_grid_seq_count = g_seq_count / pr;

  row_seq_start_idx = avg_grid_seq_count * grid_row_id;
  row_seq_end_idx = (grid_row_id == pr - 1)
                    ? g_seq_count : row_seq_start_idx + avg_grid_seq_count;
  // Subtract 1 to get the zero based inclusive end index;
  --row_seq_end_idx;

  col_seq_start_idx = avg_grid_seq_count * grid_col_id;
  col_seq_end_idx = (grid_col_id == pc - 1)
                    ? g_seq_count : col_seq_start_idx + avg_grid_seq_count;
  // Subtract 1 to get the zero based inclusive end index;
  --col_seq_end_idx;

  /*! Whose sequences will I need */

  // my_nbrs can include me as well. We'll, of course, not send
  // messages to self but I need to know what's needed from locally
  // available squences to satisfy the grid cell I am responsible for.

  find_nbrs(pr, grid_row_id, avg_l_seq_count, avg_grid_seq_count,
            row_seq_start_idx, 1, my_nbrs);
  if (grid_row_id != grid_col_id) {
    find_nbrs(pc, grid_col_id, avg_l_seq_count, avg_grid_seq_count,
              col_seq_start_idx, 0, my_nbrs);
  }

#ifndef NDEBUG
  {
    std::string title = "My neighbors";
    std::string msg;
    for (auto &nbr : my_nbrs) {
      msg += +" (rc_flag: " + std::to_string(nbr.rc_flag)
             + ", nbr_rank: " + std::to_string(nbr.nbr_rank)
             + ", nbr_seq_start_idx: " + std::to_string(nbr.nbr_seq_start_idx)
             + ", nbr_seq_end_idx: " + std::to_string(nbr.nbr_seq_end_idx) +
             ")\n";
    }
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

  // Receive neighbors are my neighbors excluding me, so receive neighbor
  // indexes will keep track of the indexes of my neighbors to receive from
  recv_nbrs_count = 0;
  for (int i = 0; i < my_nbrs.size(); ++i) {
    if (my_nbrs[i].nbr_rank == parops->world_proc_rank)
    {
      // No need to receive from self
      continue;
    }
    recv_nbrs_idxs.push_back(i);
    ++recv_nbrs_count;
  }

  /*! Issue receive requests to get the buffer lengths from my neighbors
  * for the sequences I need from them.
  */
  auto *recv_nbrs_buff_lengths_reqs = new MPI_Request[recv_nbrs_count];
  recv_nbrs_buff_lengths = new uint64_t[recv_nbrs_count];

  for (int i = 0; i < recv_nbrs_count; ++i)
  {

    MPI_Irecv(recv_nbrs_buff_lengths + i, 1, MPI_UINT64_T,
              my_nbrs[recv_nbrs_idxs[i]].nbr_rank,
              99+my_nbrs[recv_nbrs_idxs[i]].rc_flag,
              MPI_COMM_WORLD, recv_nbrs_buff_lengths_reqs + i);
  }

  /*! Who will need my sequences */
  int block_length[5] = {1, 1, 1, 1, 1};
  MPI_Aint displacement[5] = {offsetof(NbrData, rc_flag),
                              offsetof(NbrData, owner_rank),
                              offsetof(NbrData, nbr_rank),
                              offsetof(NbrData, nbr_seq_start_idx),
                              offsetof(NbrData, nbr_seq_end_idx)};
  MPI_Datatype types[] = {MPI_UNSIGNED_SHORT, MPI_INT, MPI_INT, MPI_UINT64_T,
                          MPI_UINT64_T};
  MPI_Datatype MPI_NbrData;
  MPI_Type_create_struct(5, block_length, displacement, types, &MPI_NbrData);
  MPI_Type_commit(&MPI_NbrData);

  /*! Number of neighbors will surely fit in int
   * and moreover we need this to be int to be used with
   * allgatherv, which didn't like when counts and displacements
   * were uint_64 arrays */
  int my_nbrs_count = static_cast<int>(my_nbrs.size());
  auto *all_nbrs_counts = new int[parops->world_procs_count];
  MPI_Allgather(&my_nbrs_count, 1, MPI_INT, all_nbrs_counts, 1, MPI_INT,
                MPI_COMM_WORLD);

  int all_nbrs_count = 0;
  auto *all_nbrs_displas = new int[parops->world_procs_count];
  all_nbrs_displas[0] = 0;
  for (int i = 0; i < parops->world_procs_count; ++i) {
    all_nbrs_count += all_nbrs_counts[i];
    if (i > 0) {
      all_nbrs_displas[i] = all_nbrs_displas[i - 1] + all_nbrs_counts[i - 1];
    }
  }

  auto *all_nbrs = new NbrData[all_nbrs_count];
  MPI_Allgatherv(&my_nbrs[0], my_nbrs_count, MPI_NbrData, all_nbrs,
                 all_nbrs_counts, all_nbrs_displas, MPI_NbrData,
                 MPI_COMM_WORLD);

#ifndef NDEBUG
  {
    std::string title = "All neighbor counts and displacements at the root";
    std::string msg;
    for (int i = 0; i < parops->world_procs_count; ++i) {
      msg += "(rank: " + std::to_string(i)
             + " count: " + std::to_string(all_nbrs_counts[i])
             + " displ:" + std::to_string(all_nbrs_displas[i]) + ")\n";
    }
    TraceUtils::print_msg_on_rank(title, msg, parops, 0);
  }
#endif

#ifndef NDEBUG
  {
    std::string title = "Allgathered neighbors at the root";
    std::string msg;
    int rank = 0;
    for (int i = 0; i < all_nbrs_count; ++i) {
      NbrData *n = all_nbrs + i;
      msg += "rank: " + std::to_string(rank)
             + " (rc_flag: " + std::to_string(n->rc_flag)
             + ", owner_rank: " + std::to_string(n->owner_rank)
             + ", nbr_rank: " + std::to_string(n->nbr_rank)
             + ", nbr_seq_start_idx: " + std::to_string(n->nbr_seq_start_idx)
             + ", nbr_seq_end_idx: " + std::to_string(n->nbr_seq_end_idx) +
             ")\n";
      if (i + 1 == all_nbrs_displas[rank + 1]) {
        ++rank;
      }
    }
    TraceUtils::print_msg_on_rank(title, msg, parops, 0);
  }
#endif

  std::vector<uint64_t> send_lengths;
  std::vector<uint64_t> send_start_offsets;
  std::vector<int> to_nbrs_idxs;
  to_nbrs_count = 0;
  std::sort(all_nbrs, all_nbrs + all_nbrs_count,
            [](const NbrData &a, const NbrData &b) -> bool {
              return a.nbr_rank < b.nbr_rank;
            });
  for (int i = 0; i < all_nbrs_count; ++i) 
  {
    if (all_nbrs[i].nbr_rank != parops->world_proc_rank)
    {
      // Not requesting my sequences, so continue.
      continue;
    }

    if (all_nbrs[i].owner_rank == parops->world_proc_rank) {
      // Available locally (me requesting my sequences), so continue.
      continue;
    }
    NbrData nbr = all_nbrs[i];
    uint64_t len, start_offset, end_offset;
    fd->buffer_size(nbr.nbr_seq_start_idx, nbr.nbr_seq_end_idx,
                    len, start_offset, end_offset);

    send_lengths.push_back(len);
    send_start_offsets.push_back(start_offset);
    to_nbrs_idxs.push_back(i);
    ++to_nbrs_count;
  }

  auto *to_nbrs_send_reqs = new MPI_Request[to_nbrs_count];
  for (int i = 0; i < to_nbrs_count; ++i) {
    MPI_Isend(&send_lengths[i], 1, MPI_UINT64_T,
              all_nbrs[to_nbrs_idxs[i]].owner_rank,
              99+all_nbrs[to_nbrs_idxs[i]].rc_flag, MPI_COMM_WORLD,
              to_nbrs_send_reqs + i);
  }

  /* Wait for length transfers */
  auto recv_stats = new MPI_Status[recv_nbrs_count];
  auto send_stats = new MPI_Status[to_nbrs_count];
  MPI_Waitall(recv_nbrs_count, recv_nbrs_buff_lengths_reqs, recv_stats);
  MPI_Waitall(to_nbrs_count, to_nbrs_send_reqs, send_stats);

#ifndef NDEBUG
  {
    std::string title = "My recv neighbors (excludes me) buff lengths";
    std::string msg;
    for (int i = 0; i < recv_nbrs_count; ++i) {
      NbrData nbr = my_nbrs[recv_nbrs_idxs[i]];
      msg += +" (rc_flag: " + std::to_string(nbr.rc_flag)
             + ", nbr_rank: " + std::to_string(nbr.nbr_rank)
             + ", nbr_seq_start_idx: " + std::to_string(nbr.nbr_seq_start_idx)
             + ", nbr_seq_end_idx: " + std::to_string(nbr.nbr_seq_end_idx) +
             +", nbr buff length: " +
             std::to_string(recv_nbrs_buff_lengths[i]) +
             ")\n";
    }
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

  // TODO - do multiple send/recvs below if buff length > int max

#ifndef NDEBUG
  {
    std::string title = "recv_nbrs and to_nbrs count";
    std::string msg = "(rnc: " + std::to_string(recv_nbrs_count) + " tnc: " +
                      std::to_string(to_nbrs_count) + ")";
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

  recv_nbrs_buffs = new char *[recv_nbrs_count];
  for (int i = 0; i < recv_nbrs_count; ++i) {
    recv_nbrs_buffs[i] = new char[recv_nbrs_buff_lengths[i]];
  }
  recv_nbrs_buffs_reqs = new MPI_Request[recv_nbrs_count];
  recv_nbrs_buffs_stats = new MPI_Status[recv_nbrs_count];

  for (int i = 0; i < recv_nbrs_count; ++i)
  {
    MPI_Irecv(recv_nbrs_buffs[i], static_cast<int>(recv_nbrs_buff_lengths[i]),
              MPI_CHAR, my_nbrs[recv_nbrs_idxs[i]].nbr_rank,
              77+my_nbrs[recv_nbrs_idxs[i]].rc_flag, MPI_COMM_WORLD,
              recv_nbrs_buffs_reqs + i);
  }

  to_nbrs_buffs_reqs = new MPI_Request[to_nbrs_count];
  to_nbrs_buffs_stat = new MPI_Status[to_nbrs_count];
  for (int i = 0; i < to_nbrs_count; ++i) {
    MPI_Isend(fd->buffer() + send_start_offsets[i],
              static_cast<int>(send_lengths[i]), MPI_CHAR,
              all_nbrs[to_nbrs_idxs[i]].owner_rank,
              77+all_nbrs[to_nbrs_idxs[i]].rc_flag, MPI_COMM_WORLD,
              to_nbrs_buffs_reqs + i);
  }


  delete[](to_nbrs_send_reqs);
  delete[](all_nbrs);
  delete[](all_nbrs_displas);
  delete[](all_nbrs_counts);
  MPI_Type_free(&MPI_NbrData);
  delete[](recv_nbrs_buff_lengths_reqs);

}

void DistributedFastaData::find_nbrs(const int grid_rc_procs_count,
                                     const int grid_rc_id,
                                     const uint64_t avg_l_seq_count,
                                     const uint64_t avg_grid_seq_count,
                                     const uint64_t rc_seq_start_idx,
                                     const ushort rc_flag,
                                     std::vector<NbrData> &my_nbrs)
{
  // Note, this rank might not have the sequence, if so we'll search further.
  int start_rank = static_cast<int>((rc_seq_start_idx + 1) / avg_l_seq_count);

  // GGGG: problem here when too many processes for a given file, long reads destroy the logic

  while (g_seq_offsets[start_rank] > rc_seq_start_idx)
  {
    /*! This loop is unlikely to happen unless the ranks above the original
     * start rank were heavily pruned during sequence removal step, which
     * removes sequences in lenght less than k-mer length.
     */
    --start_rank;
  }

  assert(start_rank >= 0 && start_rank < parops->world_procs_count);

  while (g_seq_offsets[start_rank] + l_seq_counts[start_rank] <
         rc_seq_start_idx)
  {
    /*! Another highly unlikely loop. This could happen if the above start_rank
     * contains much less number of sequences than the avg_l_seq_count due to
     * sequence pruning, which means we have to search in the ranks ahead of it
     * to find the sequence we are interested in. */
    ++start_rank;
  }

  assert(start_rank >= 0 && start_rank < parops->world_procs_count);

  uint64_t rc_seq_count_needed = avg_grid_seq_count;
  if (grid_rc_id == grid_rc_procs_count - 1) {
    /*! Last rank in the row or column world, so it might
     * need more than the avg_grid_seq_count */
    rc_seq_count_needed = g_seq_count - (avg_grid_seq_count * grid_rc_id);
  }

  int nbr_rank;
  uint64_t nbr_seq_start_idx, nbr_seq_end_idx;
  uint64_t count = 0;
  uint64_t seq_start_idx = rc_seq_start_idx;
  while (count < rc_seq_count_needed)
  {


    nbr_rank = start_rank;
    nbr_seq_start_idx = seq_start_idx - g_seq_offsets[start_rank];
    /* See how many sequences to grab from start_rank */
    uint64_t remaining_needed = rc_seq_count_needed - count;
    /* remaining_in_start_rank includes the current
     * seq pointed by rc_seq_start_idx */
    uint64_t remaining_in_start_rank = l_seq_counts[start_rank] -
                                       (seq_start_idx -
                                        g_seq_offsets[start_rank]);
    if (remaining_needed >= remaining_in_start_rank)
    {
      count += remaining_in_start_rank;
      nbr_seq_end_idx = ((seq_start_idx + remaining_in_start_rank - 1) -
                         g_seq_offsets[start_rank]);
    }
    else
    {
      count += remaining_needed;
      nbr_seq_end_idx =
        (seq_start_idx + remaining_needed - 1) - g_seq_offsets[start_rank];
    }

    my_nbrs.emplace_back(rc_flag, parops->world_proc_rank, nbr_rank,
                         nbr_seq_start_idx, nbr_seq_end_idx);
    ++start_rank;
    seq_start_idx = g_seq_offsets[start_rank];
  }

  assert(count == rc_seq_count_needed);
}

FastaData *DistributedFastaData::lfd() {
  return fd;
}

uint64_t DistributedFastaData::global_count() {
  return g_seq_count;
}

uint64_t DistributedFastaData::global_start_idx() {
  return g_seq_offsets[parops->world_proc_rank];
}

bool DistributedFastaData::is_ready() {
  return ready;
}

bool DistributedFastaData::is_diagonal() {
  return is_diagonal_cell;
}

void
DistributedFastaData::push_seqs(int rc_flag, FastaData *fd, uint64_t seqs_count,
                                uint64_t seq_start_idx) {
  ushort len;
  uint64_t start_offset, end_offset_inclusive;
  for (uint64_t i = 0; i < seqs_count; ++i) {
    char *buff = fd->get_sequence(seq_start_idx + i, len, start_offset,
                                  end_offset_inclusive);
    seqan::Dna5String *seq = new seqan::Dna5String(buff + start_offset, len);
    if (rc_flag == 1) {
      /*! grid row sequence */
      row_seqs.push_back(seq);
    } else {
      /*! grid col sequence */
      col_seqs.push_back(seq);
    }
  }
}

void DistributedFastaData::wait() {
  tp->times["StartDfd:MPI_Waitall(seqs)"] = std::chrono::system_clock::now();
  MPI_Waitall(recv_nbrs_count, recv_nbrs_buffs_reqs, recv_nbrs_buffs_stats);
  MPI_Waitall(to_nbrs_count, to_nbrs_buffs_reqs, to_nbrs_buffs_stat);
  tp->times["EndDfd:MPI_Waitall(seqs)"] = std::chrono::system_clock::now();

#ifndef NDEBUG
  {
    std::string title = "First character of received data from each neighbor";
    std::string msg;
    for (int i = 0; i < recv_nbrs_count; ++i) {
      msg += std::string(recv_nbrs_buffs[i], 1);
    }
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

#ifndef NDEBUG
  std::string title = "Received neighbor data";
  std::string msg;
#endif

  // tp->times["StartDfd:ExtractRecvSeqs"] = std::chrono::system_clock::now();
  int recv_nbr_idx = 0;
  for (auto &nbr : my_nbrs) {
    uint64_t nbr_seqs_count = (nbr.nbr_seq_end_idx - nbr.nbr_seq_start_idx) + 1;
    if (nbr.nbr_rank == parops->world_proc_rank) {
      /*! Local data, so create SeqAn sequences from <tt>fd</tt> */
      push_seqs(nbr.rc_flag, fd, nbr_seqs_count, nbr.nbr_seq_start_idx);
    } else {
      /*! Foreign data in a received buffer, so create a FastaData instance
       * and create SeqAn sequences from it. For foreign data, char buffers
       * only contain the required sequences, so sequence start index is 0.
       */
      uint64_t recv_nbr_l_end = recv_nbrs_buff_lengths[recv_nbr_idx] - 1;
#ifndef NDEBUG
      {
        for (int i = 0; i <= recv_nbr_l_end; ++i) {
          msg += recv_nbrs_buffs[recv_nbr_idx][i];
        }
        msg += "\n";

      }
#endif
      FastaData *recv_fd = new FastaData(recv_nbrs_buffs[recv_nbr_idx], k, 0,
                                         recv_nbr_l_end, tp, tu);

      push_seqs(nbr.rc_flag, recv_fd, nbr_seqs_count, 0);
      ++recv_nbr_idx;
    }
  }

#ifndef NDEBUG
  TraceUtils::print_msg(title, msg, parops);
#endif

  if (is_diagonal_cell) {
    /*! Diagonal cell, so col_seqs should point to the same sequences
     * as row_seqs. Also, it should not have received any col sequences
     * at this point */
    assert(col_seqs.empty());
    col_seqs.assign(row_seqs.begin(), row_seqs.end());
  }

  assert(row_seqs.size() == (row_seq_end_idx - row_seq_start_idx) + 1 &&
         col_seqs.size() == (col_seq_end_idx - col_seq_start_idx) + 1);

  ready = true;
  // tp->times["EndDfd:ExtractRecvSeqs"] = std::chrono::system_clock::now();
}

seqan::Dna5String *DistributedFastaData::row_seq(uint64_t l_row_idx) {
  return row_seqs[l_row_idx];
}

seqan::Dna5String *DistributedFastaData::col_seq(uint64_t l_col_idx) {
  return col_seqs[l_col_idx];
}
