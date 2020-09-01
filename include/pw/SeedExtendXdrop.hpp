// Created by Saliya Ekanayake on 2019-07-05.

#ifndef DIBELLA_SEEDEXTENDXDROP_HPP
#define DIBELLA_SEEDEXTENDXDROP_HPP

#include "PairwiseFunction.hpp"
#include "../AlignmentInfo.hpp"

//template <typename TSequenceValue, typename TSpec>
class SeedExtendXdrop : public PairwiseFunction{
public:

  SeedExtendXdrop(ScoringScheme scoring_scheme,
                  ushort seed_length, int xdrop, int seed_count);

  void
  PostAlignDecision(const AlignmentInfo& ai, bool& passed);

  void
  apply(uint64_t l_col_idx, uint64_t g_col_idx,
        uint64_t l_row_idx, uint64_t g_row_idx,
        seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
        ushort k,
        dibella::CommonKmers &cks, std::stringstream& ss) override;

  void
  apply_batch(seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
			        seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
			        uint64_t *lids,
			        uint64_t col_offset,
			        uint64_t row_offset,
			        PSpMat<dibella::CommonKmers>::Tuples &mattuples,
              ushort k,
			        std::ofstream &afs,
			        std::ofstream &lfs) override;

private:
  ScoringScheme scoring_scheme;
  ushort seed_length;
  int xdrop;
  int seed_count;
};

#endif //DIBELLA_SEEDEXTENDXDROP_HPP
