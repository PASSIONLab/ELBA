// Created by Saliya Ekanayake on 2019-07-10.

#include "../../include/pw/OverlapFinder.hpp"

OverlapFinder::OverlapFinder(const char *file, bool perform_alignment)
    : perform_alignment(perform_alignment), file(file){
  local_top_triangle_count = 0;
  local_nnz_count = 0;

}

void
OverlapFinder::apply(uint64_t l_col_idx, uint64_t g_col_idx, uint64_t l_row_idx,
                     uint64_t g_row_idx, seqan::Peptide *seq_h,
                     seqan::Peptide *seq_v, CommonKmers &cks, std::stringstream& ss) {

}


