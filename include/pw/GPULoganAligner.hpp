/* Created by Giulia Guidi on 4/14/2021 */

#ifndef ELBA_GPULOGAN_HPP
#define ELBA_GPULOGAN_HPP

#include "PairwiseFunction.hpp"
#include "../AlignmentInfo.hpp"

//template <typename TSequenceValue, typename TSpec>
class GPULoganAligner : public PairwiseFunction{
public:

  GPULoganAligner(ScoringScheme scoring_scheme,
                  ushort seed_length, int xdrop, int seed_count);

  void
  PostAlignDecision(const LoganAlignmentInfo& ai, bool& passed, float& ratioScoreOverlap,
        int& dir, int& dirT, int& sfx, int& sfxT, uint32_t& overlap, const bool noAlign, std::vector<int64_t>& ContainedSeqMyThread);

  void
  apply(uint64_t l_col_idx, uint64_t g_col_idx,
        uint64_t l_row_idx, uint64_t g_row_idx,
        seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
        ushort k,
        elba::CommonKmers &cks, std::stringstream& ss) override;

  void
  apply_batch(seqan::StringSet<seqan::Dna5String> &seqsh,
			        seqan::StringSet<seqan::Dna5String> &seqsv,
			        uint64_t *lids,
			        uint64_t col_offset,
			        uint64_t row_offset,
              PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
              std::ofstream &lfs,
              const bool noAlign,
              ushort k,
              uint64_t nreads,
              std::vector<int64_t>& ContainedSeqPerThread,
              float ratioScoreOverlap = 0.99,   // GGGG: Precomputed for error rate = 15% and default scoring matrix (1,-1,-1) (0.445 for CLR, 0.99 for CCS)
              int debugThr = 50) override;      // GGGG: Fixed threshold, this is convenient only for debugging

private:
  ScoringScheme scoring_scheme;
  ushort seed_length;
  int xdrop;
  int seed_count;
};

#endif //ELBA_GPULOGAN_HPP
