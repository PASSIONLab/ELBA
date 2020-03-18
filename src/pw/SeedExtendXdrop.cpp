// Created by Saliya Ekanayake on 2019-07-05.


#include "../../include/pw/SeedExtendXdrop.hpp"


//template<typename TSequenceValue, typename TSpec>
//SeedExtendXdrop<TSequenceValue, TSpec>::SeedExtendXdrop(
SeedExtendXdrop::SeedExtendXdrop(
    seqan::Blosum62 scoring_scheme,
    seqan::Blosum62 scoring_scheme_simple,
    ushort seed_length, int xdrop, int seed_count):
    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    scoring_scheme_simple(scoring_scheme_simple),
    seed_length(seed_length), xdrop(xdrop), seed_count(seed_count){

}

//template<typename TSequenceValue, typename TSpec>
//void SeedExtendXdrop<TSequenceValue, TSpec>::apply(
void SeedExtendXdrop::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Peptide *seq_h, seqan::Peptide *seq_v,
    distal::CommonKmers &cks, std::stringstream& ss) {

  AlignmentInfo ai[2];
  for (int count = 0; count < seed_count; ++count) {
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
    auto start_time = std::chrono::system_clock::now();
    extendSeed(seed, *seq_h, *seq_v, seqan::EXTEND_BOTH, scoring_scheme_simple,
               xdrop,
               seqan::GappedXDrop());
    auto end_time = std::chrono::system_clock::now();
    add_time("XA:extend_seed", (ms_t(end_time - start_time)).count());

    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), infix(*seq_h, beginPositionH(seed),
                                      endPositionH(seed)));
    assignSource(row(align, 1), infix(*seq_v, beginPositionV(seed),
                                      endPositionV(seed)));

    /*! Note. This aligns the extended seeds globally, NOT the original
     * two sequences.
     *
     * It seems kind of a waste to have to do the alignment
     * again after xdrop seed extension but that's the only
     * way to get the alignment info in SeqAn.
     * See https://seqan.readthedocs.io/en/master/Tutorial/Algorithms/SeedExtension.html
     */
    start_time = std::chrono::system_clock::now();
    globalAlignment(align, scoring_scheme);
    end_time = std::chrono::system_clock::now();
    add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());

    // Compute the statistics of the alignment.

    start_time = std::chrono::system_clock::now();
    computeAlignmentStats(ai[count].stats, align, scoring_scheme);
    end_time = std::chrono::system_clock::now();
    add_time("XA:compute_stats", (ms_t(end_time - start_time)).count());

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

  AlignmentInfo max_ai = ai[0];
  if (seed_count > 2) {
    max_ai = ai[0].stats.alignmentIdentity > ai[1].stats.alignmentIdentity
    ? ai[0] : ai[1];
  }

//  if (max_ai.stats.alignmentIdentity >= 30.0){
//    alignments.push_back(max_ai);
//  }


  double alen_minus_gapopens = (max_ai.stats.alignmentLength - max_ai.stats.numGapOpens) * 1.0;
  ss << g_col_idx << "," << g_row_idx << "," << max_ai.stats.alignmentIdentity
     << "," << max_ai.seq_h_length << "," << max_ai.seq_v_length
     << "," << max_ai.seq_h_seed_length  << "," << max_ai.seq_v_seed_length
     << "," << max_ai.stats.numGapOpens
     << "," << alen_minus_gapopens / max_ai.seq_h_length
     << "," << alen_minus_gapopens / max_ai.seq_v_length << std::endl;
}
