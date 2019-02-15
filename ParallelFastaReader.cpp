// Created by Saliya Ekanayake on 1/6/19.

#include <iostream>
#include <fstream>
#include <memory>
#include <cassert>
#include <limits>
#include <unordered_set>
#include "ParallelFastaReader.hpp"


void
readFasta(const char *file, ushort overlap, ushort k, int rank, int world_size,
          char *&buff, uint64_t &l_start, uint64_t &l_end) {
  MPI_File f;

  int err = MPI_File_open(MPI_COMM_WORLD, file,
                          MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
  if (err) {
    if (rank == 0) fprintf(stderr, "Couldn't open file %s\n", file);
    MPI_Finalize();
    exit(2);
  }

  /* The idea:
   * if sequence names and contents are roughly equal
   * then dividing the file size by the number of ranks
   * should land close to l_seq_count sequences. We'll fix
   * any mismatch later
   */

  /* Thanks, Rob Lathem for making my life a bit easy here
   * https://stackoverflow.com/a/12942718 */
  MPI_Offset g_start;
  uint64_t l_chunk_size;

  /* read in relevant chunk of file into "chunk",
   * which starts at location in the file g_start
   * and has size l_chunk_size
   */

  MPI_Offset g_end;
  MPI_Offset file_size;

  /* figure out who reads what */
  MPI_File_get_size(f, &file_size);
  file_size--;  /* get rid of text file eof */
  l_chunk_size = static_cast<uint64_t >(file_size / world_size);
  g_start = rank * l_chunk_size;
  g_end = g_start + l_chunk_size - 1;
  if (rank == world_size - 1) g_end = file_size - 1;

  /* add overlap to the end of everyone's chunk except last proc... */
  if (rank != world_size - 1)
    g_end += overlap;

  l_chunk_size = static_cast<uint64_t >(g_end - g_start + 1);


  buff = static_cast<char *>(malloc((l_chunk_size + 1) * sizeof(char)));

  /* everyone reads in their part */
  // TODO - Saliya: fix if l_chunk_size > int max. Read multiple times
  assert(l_chunk_size <= std::numeric_limits<int>::max());
  MPI_File_read_at_all(f, g_start, buff, static_cast<int>(l_chunk_size),
                       MPI_CHAR, MPI_STATUS_IGNORE);
  buff[l_chunk_size] = '\0';


  /*
   * everyone calculate what their start and end *really* are by going
   * from the first newline after start to the first newline after the
   * overlap region starts (eg, after end - overlap + 1)
   */

  l_start = 0, l_end = l_chunk_size - 1;
  if (rank != 0) {
    while (buff[l_start] != '>') l_start++;
  }
  if (rank != world_size - 1) {
    l_end -= overlap;
    while (buff[l_end] != '>') l_end++;
    // minus 2 because we don't need '>' as well as the '\n' before that
    l_end -= 2;
  }
}

std::pair<uint64_t, uint64_t>
read_local_sequences(char *buff, uint64_t l_start, uint64_t l_end, ushort k,
                     uvec_64 *id_starts, uvec_64 *seq_starts) {

  uint64_t l_seq_count = 0;
  char c;
  bool in_name = false;
  bool in_seq = false;
  ushort seq_len = 0;
  /*! No character count. This includes new line and * characters
   * It also includes entire sequences that are less than k-mer length*/
  ushort nc_count = 0;
  uint64_t idx;
  /*! Assume the FASTA content is valid */
  for (uint64_t i = l_start; i <= l_end; ++i) {
    c = buff[i];
    idx = i - nc_count;
    buff[idx] = c;
    if (c == '>') {
      id_starts->push_back(idx);
      seq_len = 0;
      in_name = true;
      in_seq = false;
      ++l_seq_count;
    } else if (c == '\n') {
      if (in_name && i + 1 <= l_end) {
        seq_starts->push_back(idx + 1);
        in_name = false;
        in_seq = true;
      } else if (in_seq && i + 1 <= l_end) {
        if (buff[i + 1] != '>') {
          ++nc_count;
        } else if (buff[i + 1] == '>') {
          if (seq_len < k) {
            ++seq_len; // capture the new line character too for removal
            uint64_t seq_id_start = id_starts->back();
            uint64_t seq_start = seq_starts->back();
            /*! + 1 to capture the new line between
             * sequence id and sequence start */
            nc_count += seq_len + (seq_start - seq_id_start) + 1;
            id_starts->pop_back();
            seq_starts->pop_back();
            --l_seq_count;
          }
        }
      }
    } else if (c == '*') {
      if (in_seq) {
        ++nc_count;
      }
    } else {
      if (in_seq) {
        ++seq_len;
      }
    }
  }

  // Remove the last sequence as well if it's shorter than k
  if (seq_len < k) {
    // The last sequence doesn't have a new line at the end
    // as our l_end points to the last real character of the block.
    uint64_t seq_id_start = id_starts->back();
    uint64_t seq_start = seq_starts->back();
    // + 1 to capture the new line between sequence id and sequence start
    nc_count += seq_len + (seq_start - seq_id_start) + 1;
    id_starts->pop_back();
    seq_starts->pop_back();
    --l_seq_count;
  }

  return {l_seq_count, nc_count};
}

std::unique_ptr<DistributedFastaData>
ParallelFastaReader
::read_and_distribute_fasta(const char *file, ushort overlap,
                            ushort k,
                            const std::shared_ptr<ParallelOps> &parops) {

  char *buff;
  uint64_t l_start, l_end;
  readFasta(file, overlap, k, parops->world_proc_rank,
            parops->world_procs_count, buff, l_start, l_end);

  auto *id_starts = new uvec_64();
  auto *seq_starts = new uvec_64();
  auto ret = read_local_sequences(buff, l_start, l_end, k,
                                  id_starts, seq_starts);

  uint64_t l_seq_count = ret.first;
  uint64_t nc_count = ret.second;

  /*! Things can go wrong if you don't end up having at least one sequence,
   * which is unlikely unless the total number of sequences are close to
   * the total number of processes.
   */
  assert(l_seq_count > 0);

  auto *l_seq_counts = new uint64_t[parops->world_procs_count];
  MPI_Allgather(&l_seq_count, 1, MPI_UINT64_T, l_seq_counts,
                1, MPI_UINT64_T, MPI_COMM_WORLD);

  auto *g_seq_offsets = new uint64_t[parops->world_procs_count];
  g_seq_offsets[0] = 0;
  for (int i = 1; i < parops->world_procs_count; ++i) {
    g_seq_offsets[i] = g_seq_offsets[i - 1] + l_seq_counts[i - 1];
  }

  /*! Total sequence count. This could be different from the original
   * input sequence count as we remove sequences less than k-mer length.
   */
  uint64_t g_seq_count = g_seq_offsets[parops->world_procs_count - 1] +
                         l_seq_counts[parops->world_procs_count - 1];

  std::unique_ptr<DistributedFastaData> dfd =
    std::make_unique<DistributedFastaData>(
      buff, l_seq_count, g_seq_count, l_start, (l_end - nc_count),
      g_seq_offsets[parops->world_proc_rank], id_starts, seq_starts);

  // TODO - Continue from here after finishing find_grid_seqs
  std::pair<char *, char *> grid_seqs
    = find_grid_seqs(g_seq_count, g_seq_offsets, l_seq_count, l_seq_counts,
                     dfd, parops);

  // TODO - DONE (see above) change this constructor to use computed values including nc_count
//  std::unique_ptr<DistributedFastaData> dfd =
//    std::make_unique<DistributedFastaData>(buff, l_start, l_end - nc_count, k);

  delete[](l_seq_counts);
  delete[](g_seq_offsets);
  return dfd;
}

std::pair<char *, char *>
ParallelFastaReader::find_grid_seqs(uint64_t g_seq_count,
                                    const uint64_t *g_seq_offsets,
                                    uint64_t l_seq_count,
                                    const uint64_t *l_seq_counts,
                                    std::unique_ptr<DistributedFastaData> &dfd,
                                    const std::shared_ptr<ParallelOps> &parops) {
  uint64_t avg_l_seq_count = g_seq_count / parops->world_procs_count;
  uint64_t pr, pc;
  pr = pc = static_cast<uint64_t>(sqrt(parops->world_procs_count));
  uint64_t avg_grid_seq_count = g_seq_count / pr;

  /*! It might seem wrong to use rank_in_col for row_id and rank_in_row for
   * col_id but it is correct. These are ranks when the column and row worlds
   * are viewed separately.
   */
  auto grid_row_id = static_cast<uint64_t>(parops->grid->GetRankInProcCol());
  auto grid_col_id = static_cast<uint64_t>(parops->grid->GetRankInProcRow());

  uint64_t row_seq_start_idx = avg_grid_seq_count * grid_row_id;
  uint64_t col_seq_start_idx = avg_grid_seq_count * grid_col_id;

  /*! Whose sequences will I need */

  // my_nbrs can include me as well. We'll, of course, not send
  // messages to self but I need to know what's needed from locally
  // available squences to satisfy the grid cell I am responsible for.
  std::vector<NbrData> my_nbrs;
  find_nbrs(pr, grid_row_id, g_seq_count, avg_l_seq_count, avg_grid_seq_count,
            row_seq_start_idx, 1, g_seq_offsets, l_seq_counts, parops, my_nbrs);
  if (grid_row_id != grid_col_id) {
    find_nbrs(pc, grid_col_id, g_seq_count, avg_l_seq_count, avg_grid_seq_count,
              col_seq_start_idx, 0, g_seq_offsets, l_seq_counts, parops,
              my_nbrs);
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
    DebugUtils::print_msg(title, msg, parops);
  }
#endif


  // Receive neighbors are my neighbors excluding me, so receive neighbor
  // indexes will keep track of the indexes of my neighbors to receive from
  std::vector<int> recv_nbrs_idxs;
  int recv_nbrs_count = 0;
  for (int i = 0; i < my_nbrs.size(); ++i) {
    if (my_nbrs[i].nbr_rank == parops->world_proc_rank) {
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
  auto recv_nbrs_buff_lengths = new uint64_t[recv_nbrs_count];
  for (int i = 0; i < recv_nbrs_count; ++i) {
    MPI_Irecv(recv_nbrs_buff_lengths + i, 1, MPI_UINT64_T,
              my_nbrs[recv_nbrs_idxs[i]].nbr_rank, 99,
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
    DebugUtils::print_msg_on_rank(title, msg, parops, 0);
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
    DebugUtils::print_msg_on_rank(title, msg, parops, 0);
  }
#endif

  std::vector<uint64_t> send_lengths;
  std::vector<uint64_t> send_start_offsets;
  std::vector<int> to_nbrs_idxs;
  int to_nbrs_count = 0;
  std::sort(all_nbrs, all_nbrs + all_nbrs_count,
            [](const NbrData &a, const NbrData &b) -> bool {
              return a.nbr_rank < b.nbr_rank;
            });
  for (int i = 0; i < all_nbrs_count; ++i) {
    if (all_nbrs[i].nbr_rank != parops->world_proc_rank) {
      // Not requesting my sequences, so continue.
      continue;
    }

    if (all_nbrs[i].owner_rank == parops->world_proc_rank) {
      // Available locally (me requesting my sequences), so continue.
      continue;
    }
    NbrData nbr = all_nbrs[i];
    uint64_t len, start_offset, end_offset;
    dfd->buffer_size(nbr.nbr_seq_start_idx, nbr.nbr_seq_end_idx,
                     len, start_offset, end_offset);

    send_lengths.push_back(len);
    send_start_offsets.push_back(start_offset);
    to_nbrs_idxs.push_back(i);
    ++to_nbrs_count;
  }

  auto *to_nbrs_send_reqs = new MPI_Request[to_nbrs_count];
  for (int i = 0; i < to_nbrs_count; ++i) {
    MPI_Isend(&send_lengths[i], 1, MPI_UINT64_T,
              all_nbrs[to_nbrs_idxs[i]].owner_rank, 99, MPI_COMM_WORLD,
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
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

  // TODO - do multiple send/recvs if buff length > int max

  /* Now, I know the recv buffer lengths from my recv neighbors, so let's
   * do the actual exchange */
//  char **recv_nbrs_data_buffs = new char *[recv_nbrs_count];
//  for (int i = 0; i < recv_nbrs_count; ++i){
//    recv_nbrs_data_buffs[i] = new char[static_cast<int>(recv_nbrs_buff_lengths[i])];
//  }
//  auto *recv_nbrs_buff_reqs = new MPI_Request[recv_nbrs_count];
//  auto *recv_nbrs_buff_stats = new MPI_Status[recv_nbrs_count];
//
//  for (int i = 0; i < recv_nbrs_count; ++i) {
////    MPI_Irecv(recv_nbrs_data_buffs[i],
////              static_cast<int>(recv_nbrs_buff_lengths[i]), MPI_CHAR,
////              my_nbrs[recv_nbrs_idxs[i]].nbr_rank, 98, MPI_COMM_WORLD,
////              recv_nbrs_buff_reqs + i);
//    MPI_Irecv(recv_nbrs_data_buffs[i],
//              1, MPI_CHAR,
//              my_nbrs[recv_nbrs_idxs[i]].nbr_rank, 98, MPI_COMM_WORLD,
//              recv_nbrs_buff_reqs + i);
//  }
//
//  auto *to_nbrs_buff_reqs = new MPI_Request[recv_nbrs_count];
//  auto *to_nbrs_buff_stats = new MPI_Status[recv_nbrs_count];
//  for (int i = 0; i < to_nbrs_count; ++i) {
////    MPI_Isend(dfd->buffer() + send_start_offsets[i],
////              static_cast<int>(send_lengths[i]), MPI_CHAR,
////              all_nbrs[to_nbrs_idxs[i]].owner_rank, 98, MPI_COMM_WORLD,
////              to_nbrs_buff_reqs+i);
//    char c = static_cast<char>(parops->world_proc_rank + 65);
//    MPI_Isend(&c,
//              1, MPI_CHAR,
//              all_nbrs[to_nbrs_idxs[i]].owner_rank, 98, MPI_COMM_WORLD,
//              to_nbrs_buff_reqs+i);
//  }


#ifndef NDEBUG
  {
    std::string title = "recv_nbrs and to_nbrs count";
    std::string msg = "(rnc: " + std::to_string(recv_nbrs_count) + " tnc: " + std::to_string(to_nbrs_count) + ")";
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

// TODO - Debug
  char *debug_recvs = new char[recv_nbrs_count];
  auto *debug_recv_reqs = new MPI_Request[recv_nbrs_count];
  auto *debug_recv_stats = new MPI_Status[recv_nbrs_count];
  for (int i = 0; i < recv_nbrs_count; ++i){
    MPI_Irecv(debug_recvs+i, 1, MPI_CHAR, my_nbrs[recv_nbrs_idxs[i]].nbr_rank, 77, MPI_COMM_WORLD, debug_recv_reqs+i);
  }

  char debug_v = static_cast<char>(parops->world_proc_rank + 48);
  auto *debug_send_reqs = new MPI_Request[to_nbrs_count];
  auto *debug_send_stats = new MPI_Status[to_nbrs_count];
  for (int i = 0; i < to_nbrs_count; ++i){
    MPI_Isend(&debug_v, 1, MPI_CHAR, all_nbrs[to_nbrs_idxs[i]].owner_rank, 77, MPI_COMM_WORLD, debug_send_reqs+i);
  }

  MPI_Waitall(recv_nbrs_count, debug_recv_reqs, debug_recv_stats);
  MPI_Waitall(to_nbrs_count, debug_send_reqs, debug_send_stats);

#ifndef NDEBUG
  {
    std::string title = "debug";
    std::string msg(debug_recvs, recv_nbrs_count);
    DebugUtils::print_msg(title, msg, parops);
  }
#endif

  delete[](debug_send_stats);
  delete[](debug_send_reqs);
  delete[](debug_recv_stats);
  delete[](debug_recv_reqs);
  delete[](debug_recvs);

//  MPI_Waitall(to_nbrs_count, to_nbrs_buff_reqs, to_nbrs_buff_stats);
//  MPI_Waitall(recv_nbrs_count, recv_nbrs_buff_reqs, recv_nbrs_buff_stats);
//
//  delete[](to_nbrs_buff_stats);
//  delete[](to_nbrs_buff_reqs);
//  delete[](recv_nbrs_buff_stats);
//  delete[](recv_nbrs_buff_reqs);
//  for (int i = 0; i < recv_nbrs_count; ++i) {
//    delete[](recv_nbrs_data_buffs[i]);
//  }
//  delete[](recv_nbrs_data_buffs);
  delete[](to_nbrs_send_reqs);
  delete[](all_nbrs);
  delete[](all_nbrs_displas);
  delete[](all_nbrs_counts);
  delete[](recv_nbrs_buff_lengths_reqs);
  delete[](recv_nbrs_buff_lengths);
  MPI_Type_free(&MPI_NbrData);

  return std::pair<char *, char *>();
}

void
ParallelFastaReader::find_nbrs(
  const uint64_t grid_rc_procs_count,
  const uint64_t grid_rc_id,
  const uint64_t g_seq_count,
  const uint64_t avg_l_seq_count,
  const uint64_t avg_grid_seq_count,
  const uint64_t rc_seq_start_idx,
  const ushort rc_flag,
  const uint64_t *g_seq_offsets,
  const uint64_t *l_seq_counts,
  const std::shared_ptr<ParallelOps> &parops,
  std::vector<NbrData> &nbrs
) {
  // Note, this rank might not have the sequence, if so we'll search further.
  int start_rank
    = static_cast<int>((rc_seq_start_idx + 1) / avg_l_seq_count);
  while (g_seq_offsets[start_rank] > rc_seq_start_idx) {
    /*! This loop is unlikely to happen unless the ranks above the original
     * start rank were heavily pruned during sequence removal step, which
     * removes sequences in lenght less than k-mer length.
     */
    --start_rank;
  }

  assert(start_rank >= 0 && start_rank < parops->world_procs_count);

  while (g_seq_offsets[start_rank] + l_seq_counts[start_rank] <
         rc_seq_start_idx) {
    /* Another highly unlikely loop. This could happen if the above start_rank
     * contains much less number of sequences than the avg_l_seq_count due to
     * sequence pruning, which means we have to search in the ranks ahead of it
     * to find the sequence we are interested in. */
    ++start_rank;
  }

  /*{
    *//* TODO - debug *//*
    std::string str = "(start_rank:" + std::to_string(start_rank)
                      + " avg_l_seq_count:" + std::to_string(avg_l_seq_count)
                      + " rc_seq_start_idx:" +
                      std::to_string(rc_seq_start_idx)
                      + ")\n";
    int flag;
    if (parops->world_proc_rank > 0)
      MPI_Recv(&flag, 1, MPI_INT, parops->world_proc_rank - 1,
               99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout << "\nRank: " << parops->world_proc_rank
              << "\n------------------------\n";
    std::cout << str;
    if (parops->world_proc_rank < parops->world_procs_count - 1) {
      MPI_Send(&flag, 1, MPI_INT, parops->world_proc_rank + 1, 99,
               MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    *//* End debug *//*
  }*/

  assert(start_rank >= 0 && start_rank < parops->world_procs_count);

  uint64_t rc_seq_count_needed = avg_grid_seq_count;
  if (grid_rc_id == grid_rc_procs_count - 1) {
    /*! Last rank in the column world, so it might
     * need more than the avg_grid_seq_count */
    rc_seq_count_needed = g_seq_count - (avg_grid_seq_count * grid_rc_id);
  }

  int nbr_rank;
  uint64_t nbr_seq_start_idx, nbr_seq_end_idx;
  uint64_t count = 0;
  uint64_t seq_start_idx = rc_seq_start_idx;
  while (count < rc_seq_count_needed) {
    nbr_rank = start_rank;
    nbr_seq_start_idx = seq_start_idx - g_seq_offsets[start_rank];
    /* See how many sequences to grab from start_rank */
    uint64_t remaining_needed = rc_seq_count_needed - count;
    /* remaining_in_start_rank includes the current
     * seq pointed by rc_seq_start_idx */
    uint64_t remaining_in_start_rank = l_seq_counts[start_rank] -
                                       (seq_start_idx -
                                        g_seq_offsets[start_rank]);
    if (remaining_needed >= remaining_in_start_rank) {
      count += remaining_in_start_rank;
      nbr_seq_end_idx = ((seq_start_idx + remaining_in_start_rank - 1) -
                         g_seq_offsets[start_rank]);
    } else {
      count += remaining_needed;
      nbr_seq_end_idx =
        (seq_start_idx + remaining_needed - 1) - g_seq_offsets[start_rank];
    }
    nbrs.emplace_back(rc_flag, parops->world_proc_rank, nbr_rank,
                      nbr_seq_start_idx, nbr_seq_end_idx);
    ++start_rank;
    seq_start_idx = g_seq_offsets[start_rank];
  }

  assert(count == rc_seq_count_needed);
}






