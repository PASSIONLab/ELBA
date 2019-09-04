// Created by Saliya Ekanayake on 2019-07-05.

#ifndef LBL_PISA_PAIRWISEFUNCTION_HPP
#define LBL_PISA_PAIRWISEFUNCTION_HPP

#include <seqan/score.h>
#include "../Kmer.hpp"
#include "../AlignmentInfo.hpp"

class PairwiseFunction {
public:
  PairwiseFunction();
  virtual ~PairwiseFunction();

  virtual void apply(uint64_t l_col_idx, uint64_t g_col_idx,
      uint64_t l_row_idx, uint64_t g_row_idx,
      seqan::Peptide *seq_h, seqan::Peptide *seq_v,
      CommonKmers &cks) = 0;

};

#endif //LBL_PISA_PAIRWISEFUNCTION_HPP
