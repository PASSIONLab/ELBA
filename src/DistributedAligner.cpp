//
// Created by Saliya Ekanayake on 2019-02-19.
//

#include "../include/DistributedAligner.hpp"

DistributedAligner::DistributedAligner(ushort seed_length, int xdrop,
                                       int gap_open, int gap_ext,
                                       const std::shared_ptr<DistributedFastaData> dfd,
                                       PSpMat<CommonKmers>::MPI_DCCols mat,
                                       const std::shared_ptr<ParallelOps> &parops)
  : seed_length(seed_length), xdrop(xdrop), gap_open(gap_open),
    gap_ext(gap_ext),
    dfd(dfd), mat(mat), parops(parops) {

}

uint64_t DistributedAligner::align_seqs() {
  /*! There are two types of rows and columns below.
   * The sequences are arranged as an NxN matrix in
   * mat (this is not how it's stored internally).
   * This NxN is distributed over a grid of
   * sqrt(P) x sqrt (P), where P is the total number
   * of processes. Anything to do with the grid will
   * be prefixed by gr_*/

  // rows and cols in the result
  uint64_t n_rows, n_cols;
  n_rows = n_cols = dfd->global_count();
  int gr_rows = parops->grid->GetGridRows();
  int gr_cols = parops->grid->GetGridCols();

  int gr_col_idx = parops->grid->GetRankInProcRow();
  int gr_row_idx = parops->grid->GetRankInProcCol();
  uint64_t avg_rows_in_grid = n_rows / gr_rows;
  uint64_t avg_cols_in_grid = n_cols / gr_cols;
  // local submatrix
  PSpMat<CommonKmers>::DCCols *spSeq = mat.seqptr();
  // first row in this process
  uint64_t row_offset = gr_row_idx * avg_rows_in_grid;
  // first col in this process
  uint64_t col_offset = gr_col_idx * avg_cols_in_grid;

  seqan::Blosum62 blosum62(gap_ext, gap_open);
  // TODO: SeqAn can't work with affine gaps for seed extension
  seqan::Blosum62 blosum62_simple(gap_open, gap_open);

  for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit) {
    // iterate over columns
    auto l_col_idx = colit.colid(); // local numbering
    uint64_t g_col_idx = l_col_idx + col_offset;

    seqan::Peptide *seq_h = dfd->col_seq(l_col_idx);
    for (auto nzit = spSeq->begnz(colit); nzit < spSeq->endnz(colit); ++nzit) {
      auto l_row_idx = nzit.rowid();
      uint64_t g_row_idx = l_row_idx + row_offset;

      seqan::Peptide *seq_v = dfd->row_seq(l_row_idx);

      /*!
       * Note. the cells means the process grid cells.
       * We only want to compute the top triangles of any grid cell.
       * Further, we want cell diagonals only for cells that are on the
       * top half of the grid excluding the grid's main diagonal cells
       */
      if (l_col_idx < l_row_idx){
        continue;
      }

      if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx){
        continue;
      }

      CommonKmers cks = nzit.value();

      AlignmentInfo ai[2];
      for (int count = 0; count < 2; ++count) {
        // row sequence is the same thing as vertical sequence
        ushort l_row_seed_start_offset = (count == 0) ? cks.first.first
                                                      : cks.second.first;
        // col sequence is the same thing as horizontal sequence
        ushort l_col_seed_start_offset = (count == 0) ? cks.first.second
                                                      : cks.second.second;

        // Seed creation params are:
        // horizontal seed start offset, vertical seed start offset, length
        TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset,
                   seed_length);
        extendSeed(seed, *seq_h, *seq_v, seqan::EXTEND_BOTH, blosum62_simple,
                   xdrop,
                   seqan::GappedXDrop());

        seqan::Align<seqan::Peptide> align;
        resize(rows(align), 2);
        assignSource(row(align, 0), infix(*seq_h, beginPositionH(seed),
                                          endPositionH(seed)));
        assignSource(row(align, 1), infix(*seq_v, beginPositionV(seed),
                                          endPositionV(seed)));


        globalAlignment(align, blosum62);

        // Compute the statistics of the alignment.

        computeAlignmentStats(ai[count].stats, align, blosum62);
        ai[count].seq_h_length = length(*seq_h);
        ai[count].seq_v_length = length(*seq_v);
        ai[count].seq_h_seed_length = static_cast<ushort>(seed._endPositionH -
                                                          seed._beginPositionH);
        ai[count].seq_v_seed_length = static_cast<ushort>(seed._endPositionV -
                                                          seed._beginPositionV);
        ai[count].seq_h_g_idx = g_col_idx;
        ai[count].seq_v_g_idx = gr_row_idx;
      }

      /* Hard coding quality constraints for now */

      // TODO - Saliya
      // For now only keeps the largest alignment > 30% identity.
      // Incorporate length coverage restrictions later.
      AlignmentInfo max_ai =
        ai[0].stats.alignmentIdentity > ai[1].stats.alignmentIdentity
        ? ai[0] : ai[1];

//      if (max_ai.stats.alignmentIdentity < 30.0){
//        continue;
//      }
      alignments.push_back(max_ai);
    }
  }

  uint64_t local_alignments = alignments.size();
  uint64_t total_alignments = 0;
  MPI_Reduce(&local_alignments, &total_alignments, 1, MPI_UINT64_T, MPI_SUM, 0,
             MPI_COMM_WORLD);
  return total_alignments;
}

void DistributedAligner::write_overlaps(const char *file) {
  // rows and cols in the result
  uint64_t n_rows, n_cols;
  n_rows = n_cols = dfd->global_count();
  int gr_rows = parops->grid->GetGridRows();
  int gr_cols = parops->grid->GetGridCols();

  int gr_col_idx = parops->grid->GetRankInProcRow();
  int gr_row_idx = parops->grid->GetRankInProcCol();
  uint64_t avg_rows_in_grid = n_rows / gr_rows;
  uint64_t avg_cols_in_grid = n_cols / gr_cols;
  // local submatrix
  PSpMat<CommonKmers>::DCCols *spSeq = mat.seqptr();
  // first row in this process
  uint64_t row_offset = gr_row_idx * avg_rows_in_grid;
  // first col in this process
  uint64_t col_offset = gr_col_idx * avg_cols_in_grid;

  uint64_t local_nnz_count = 0;
  uint64_t local_top_triangle_count = 0;
  std::stringstream ss;
  for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit) {
    // iterate over columns
    auto l_col_idx = colit.colid(); // local numbering
    uint64_t g_col_idx = l_col_idx + col_offset;

    for (auto nzit = spSeq->begnz(colit); nzit < spSeq->endnz(colit); ++nzit) {
      auto l_row_idx = nzit.rowid();
      uint64_t g_row_idx = l_row_idx + row_offset;

      ++local_nnz_count;

      /*!
       * Note. the cells means the process grid cells.
       * We only want to compute the top triangles of any grid cell.
       * Further, we want cell diagonals only for cells that are on the
       * top half of the grid excluding the grid's main diagonal cells
       */
      if (l_col_idx < l_row_idx){
        continue;
      }

      if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx){
        continue;
      }

      ++local_top_triangle_count;
      ss << g_col_idx << "," << g_row_idx << "\n";
    }
  }

  std::string overlaps_str = ss.str();
  parops->write_file_in_parallel(file, overlaps_str);
}
