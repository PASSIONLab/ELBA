// Created by Saliya Ekanayake on 2019-07-10.

#ifndef ELBA_OVERLAPFINDER_HPP
#define ELBA_OVERLAPFINDER_HPP

#include "PairwiseFunction.hpp"

class OverlapFinder : PairwiseFunction{
public:
  OverlapFinder(const char *file, bool perform_alignment);

  void apply(uint64_t l_col_idx, uint64_t g_col_idx, uint64_t l_row_idx,
             uint64_t g_row_idx, seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
             ushort k,
             dibella::CommonKmers &cks, std::stringstream& ss) override;

private:
  bool perform_alignment;
  const char *file;

  uint64_t local_top_triangle_count;
  uint64_t local_nnz_count;
};

#endif //ELBA_OVERLAPFINDER_HPP
