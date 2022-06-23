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
  PostAlignDecision(const AlignmentInfo& ai, bool& passed, float& ratioScoreOveralap,
          int& dir, int& dirT, int& sfx, int& sfxT, uint32_t& overlap, const bool noAlign, std::vector<int64_t>& ContainedSeqPerProc);

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
              PSpMat<dibella::CommonKmers>::ref_tuples *mattuples,
              std::ofstream &lfs,
              const bool noAlign,
              ushort k,
              uint64_t nreads,
              std::vector<int64_t>& ContainedSeqPerProc,
              float ratioScoreOverlap = 0.99,    // GGGG: Precomputed for error rate = 15% and default scoring matrix (1,-1,-1) (0.445 for CLR, 0.99 for CCS)
              //float ratioScoreOverlap = 0.445,    // GGGG: Precomputed for error rate = 15% and default scoring matrix (1,-1,-1) (0.445 for CLR, 0.99 for CCS)
              int debugThr = 50) override;      // GGGG: Fixed threshold, this is convenient only for debugging

private:
  ScoringScheme scoring_scheme;
  ushort seed_length;
  int xdrop;
  int seed_count;
};

#endif //DIBELLA_SEEDEXTENDXDROP_HPP
