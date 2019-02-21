//
// Created by Saliya Ekanayake on 2019-02-19.
//

#include "DistributedAligner.hpp"

DistributedAligner::DistributedAligner(ushort seed_length, int xdrop,
                                       int gap_open, int gap_ext,
                                       const std::shared_ptr<DistributedFastaData> dfd,
                                       PSpMat<CommonKmers>::MPI_DCCols mat,
                                       const std::shared_ptr<ParallelOps> &parops)
  : seed_length(seed_length), xdrop(xdrop), gap_open(gap_open),
    gap_ext(gap_ext),
    dfd(dfd), mat(mat), parops(parops) {

}

void DistributedAligner::align() {
  /*! There are two types of rows and columns below.
   * The sequences are arranged as an NxN matrix in
   * mat (this is not how it's store internally).
   * This N is distributed over a grid of
   * sqrt(P) x sqrt (P), where P is the toal number
   * of processes. Anything to do with the grid will
   * be prefixed by gr_*/

  // rows and cols in the result
  uint64_t rows, cols;
  rows = cols = dfd->global_count();
  int gr_rows = parops->grid->GetGridRows();
  int gr_cols = parops->grid->GetGridCols();

  int gr_col_idx = parops->grid->GetRankInProcRow();
  int gr_row_idx = parops->grid->GetRankInProcCol();
  uint64_t avg_rows_in_grid = rows / gr_rows;
  uint64_t avg_cols_in_grid = cols / gr_cols;
  // local submatrix
  PSpMat<CommonKmers>::DCCols *spSeq = mat.seqptr();
  // first row in this process
  uint64_t row_offset = gr_row_idx * avg_rows_in_grid;
  // first col in this process
  uint64_t col_offset = gr_col_idx * avg_cols_in_grid;

  for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit) {
    // iterate over columns
    auto l_col_idx = colit.colid(); // local numbering
    uint64_t g_col_idx = l_col_idx + col_offset;

    seqan::Peptide seq_h = dfd->col_seq(l_col_idx);
    for (auto nzit = spSeq->begnz(colit); nzit < spSeq->endnz(colit); ++nzit) {
      auto l_row_idx = nzit.rowid();
      uint64_t g_row = l_row_idx + row_offset;

      seqan::Peptide seq_v = dfd->row_seq(l_row_idx);

      /*! Align only the top part of this cell. Further, if this is a
       * diagonal cell then avoid aligning the diagonal pairs
       */

      if (dfd->is_diagonal() && l_col_idx <= l_row_idx) {
        continue;
      }

      if (l_col_idx < l_row_idx) {
        continue;
      }

      CommonKmers cks = nzit.value();
      // row sequence is the same thing as vertical sequence
      ushort l_row_seed_start_offset = cks.first.first;
      // col sequence is the same thing as horizontal sequence
      ushort l_col_seed_start_offset = cks.first.second;

      // Seed creation params are:
      // horizontal seed start offset, vertical seed start offset, length
//      TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset, seed_length);
//      extendSeed(seed, seq_h, seq_v, seqan::EXTEND_BOTH, xdrop,
//                 seqan::GappedXDrop());







//      rf << g_row << ",";
//      cf << g_col_idx << ",";
//      vf << nzit.value().count << ",";
//      std::cout<<"r:"<<l_row_idx<<" c:"<<l_col_idx<<" v:"<<nzit.value()<<std::endl;
    }
  }

}
