// Created by Saliya Ekanayake on 2019-09-03.

#include "../../include/pw/FullAligner.hpp"

FullAligner::FullAligner(ScoringScheme scoring_scheme) :
    PairwiseFunction(),
    scoring_scheme(scoring_scheme)
	{}

void FullAligner::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
	ushort k,
    elba::CommonKmers &cks, std::stringstream &ss) {

	// ...
}

void
FullAligner::apply_batch
(
    seqan::StringSet<seqan::Dna5String> &seqsh,
	seqan::StringSet<seqan::Dna5String> &seqsv,
	uint64_t *lids,
	uint64_t col_offset,
	uint64_t row_offset,
    PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
    std::ofstream &lfs,
	const bool noAlign,
	ushort k,
	uint64_t nreads,
	std::vector<int64_t>& ContainedSeqPerProc,
    float ratioScoreOverlap,
	int debugThr
)
{
	// ...
}
