// Created by Saliya Ekanayake on 2019-07-05.

#ifndef DISTAL_SEEDEXTENDXDROP_HPP
#define DISTAL_SEEDEXTENDXDROP_HPP

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
             distal::CommonKmers &cks, std::stringstream& ss) override;

  void
  apply_batch (seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsh,
			   seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsv,
			   uint64_t *lids,
			   uint64_t col_offset,
			   uint64_t row_offset,
			   PSpMat<distal::CommonKmers>::Tuples &mattuples,
			   std::ofstream &afs,
			   std::ofstream &lfs) override;
  
  // std::vector<AlignmentInfo> alignments;

private:
//  seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme;
//  seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme_simple;
  seqan::Blosum62 scoring_scheme;
  seqan::Blosum62 scoring_scheme_simple;
  ushort seed_length;
  int xdrop;
  int seed_count;

};

#endif //DISTAL_SEEDEXTENDXDROP_HPP
