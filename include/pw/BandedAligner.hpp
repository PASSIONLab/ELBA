// Created by Saliya Ekanayake on 2019-09-03.

#ifndef DISTAL_BANDEDALIGNER_HPP
#define DISTAL_BANDEDALIGNER_HPP

#include "PairwiseFunction.hpp"
#include "../AlignmentInfo.hpp"

class BandedAligner  : public PairwiseFunction{
public:
//  SeedExtendXdrop(seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme,
//      seqan::Score<int, seqan::ScoreMatrix<TSequenceValue, TSpec>> scoring_scheme_simple,
//  ushort seed_length, int xdrop);

  BandedAligner(seqan::Blosum62 scoring_scheme, int banded_half_width);


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
  int banded_half_width;

};


#endif //DISTAL_BANDEDALIGNER_HPP
