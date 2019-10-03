// Created by Saliya Ekanayake on 2019-07-05.

#ifndef LBL_PISA_SEEDEXTENDXDROP_HPP
#define LBL_PISA_SEEDEXTENDXDROP_HPP

#include "PairwiseFunction.hpp"
#include "../AlignmentInfo.hpp"

//template <typename TSequenceValue, typename TSpec>
class SeedExtendXdrop : public PairwiseFunction{
public:
//  SeedExtendXdrop(seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme,
//      seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme_simple,
//  ushort seed_length, int xdrop);

  SeedExtendXdrop(seqan::Blosum62 scoring_scheme,
                  seqan::Blosum62 scoring_scheme_simple,
                  ushort seed_length, int xdrop, int seed_count);


  void apply(uint64_t l_col_idx, uint64_t g_col_idx,
             uint64_t l_row_idx, uint64_t g_row_idx,
             seqan::Peptide *seq_h, seqan::Peptide *seq_v,
             CommonKmers &cks, std::stringstream& ss) override;

  std::vector<AlignmentInfo> alignments;

private:
//  seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme;
//  seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme_simple;
  seqan::Blosum62 scoring_scheme;
  seqan::Blosum62 scoring_scheme_simple;
  ushort seed_length;
  int xdrop;
  int seed_count;

};

#endif //LBL_PISA_SEEDEXTENDXDROP_HPP
