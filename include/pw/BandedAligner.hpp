// Created by Saliya Ekanayake on 2019-09-03.

#ifndef DIBELLA_BANDEDALIGNER_HPP
#define DIBELLA_BANDEDALIGNER_HPP

#include "PairwiseFunction.hpp"
#include "../AlignmentInfo.hpp"

class BandedAligner : public PairwiseFunction{
public:

  BandedAligner(ScoringScheme scoring_scheme, int banded_half_width);

  void apply(uint64_t l_col_idx, uint64_t g_col_idx,
             uint64_t l_row_idx, uint64_t g_row_idx,
             seqan::Peptide *seq_h, seqan::Peptide *seq_v,
             dibella::CommonKmers &cks, std::stringstream& ss) override;

  void
  apply_batch(seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsh,
			        seqan::StringSet<seqan::Gaps<seqan::Peptide>> &seqsv,
			        uint64_t *lids,
			        uint64_t col_offset,
			        uint64_t row_offset,
			        PSpMat<dibella::CommonKmers>::Tuples &mattuples,
			        std::ofstream &afs,
			        std::ofstream &lfs) override;

private:
  ScoringScheme scoring_scheme;
  int banded_half_width;
};

#endif //DIBELLA_BANDEDALIGNER_HPP
