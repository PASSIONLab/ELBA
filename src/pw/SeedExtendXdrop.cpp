// Created by Saliya Ekanayake on 2019-07-05.


#include "../../include/pw/SeedExtendXdrop.hpp"


//template<typename TSequenceValue, typename TSpec>
//SeedExtendXdrop<TSequenceValue, TSpec>::SeedExtendXdrop(
SeedExtendXdrop::SeedExtendXdrop(
    seqan::Blosum62 scoring_scheme,
    seqan::Blosum62 scoring_scheme_simple,
    ushort seed_length, int xdrop):
    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    scoring_scheme_simple(scoring_scheme_simple),
    seed_length(seed_length), xdrop(xdrop){

}

//template<typename TSequenceValue, typename TSpec>
//void SeedExtendXdrop<TSequenceValue, TSpec>::apply(
void SeedExtendXdrop::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Peptide *seq_h, seqan::Peptide *seq_v,
    CommonKmers &cks) {

  AlignmentInfo ai[2];
  for (int count = 0; count < 2; ++count) {
    // row sequence is the same thing as vertical sequence
    ushort l_row_seed_start_offset = (count == 0) ? cks.first.first
                                                  : cks.second.first;
    // col sequence is the same thing as horizontal sequence
    ushort l_col_seed_start_offset = (count == 0) ? cks.first.second
                                                  : cks.second.second;

    // Seed creation params are:
    // horizontal seed start offset, vertical seed start offset, length
    TSeed seed(l_col_seed_start_offset, l_row_seed_start_offset,
               seed_length);
    extendSeed(seed, *seq_h, *seq_v, seqan::EXTEND_BOTH, scoring_scheme_simple,
               xdrop,
               seqan::GappedXDrop());

    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(*seq_h, beginPositionH(seed),
                                      endPositionH(seed)));
    assignSource(row(align, 1), infix(*seq_v, beginPositionV(seed),
                                      endPositionV(seed)));

    /*! Note. This aligns the extended seeds globally, NOT the original
     * two sequences.
     */
    globalAlignment(align, scoring_scheme);

    // Compute the statistics of the alignment.

    computeAlignmentStats(ai[count].stats, align, scoring_scheme);
    ai[count].seq_h_length = length(*seq_h);
    ai[count].seq_v_length = length(*seq_v);
    ai[count].seq_h_seed_length = static_cast<ushort>(seed._endPositionH -
                                                      seed._beginPositionH);
    ai[count].seq_v_seed_length = static_cast<ushort>(seed._endPositionV -
                                                      seed._beginPositionV);
    ai[count].seq_h_g_idx = g_col_idx;
    ai[count].seq_v_g_idx = g_row_idx;
  }

  /* Hard coding quality constraints for now */

  // TODO - Saliya
  // For now only keeps the largest alignment > 30% identity.
  // Incorporate length coverage restrictions later.
  AlignmentInfo max_ai =
      ai[0].stats.alignmentIdentity > ai[1].stats.alignmentIdentity
      ? ai[0] : ai[1];

//  if (max_ai.stats.alignmentIdentity >= 30.0){
    alignments.push_back(max_ai);
//  }
}
