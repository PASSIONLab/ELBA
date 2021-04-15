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
    dibella::CommonKmers &cks, std::stringstream &ss) {

	// ...
}

void
FullAligner::apply_batch
(
    seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
	seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
	uint64_t *lids,
	uint64_t col_offset,
	uint64_t row_offset,
    PSpMat<dibella::CommonKmers>::ref_tuples *mattuples,
    std::ofstream &lfs,
	const bool noAlign,
	ushort k,
	uint64_t nreads,
    float ratioScoreOverlap,
	int debugThr
)
{
	seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	uint64_t npairs = seqan::length(seqsh);
	setNumThreads(exec_policy, numThreads);

	lfs << "processing batch of size " << npairs << " with "
		<< numThreads << " threads " << std::endl;

	auto start_time = std::chrono::system_clock::now();
	
	// alignment
	localAlignment(exec_policy, seqsh, seqsv, scoring_scheme);
	
	auto end_time = std::chrono::system_clock::now();
  	add_time("FA:local_alignment", (ms_t(end_time - start_time)).count());

	start_time = std::chrono::system_clock::now();
	
	// stats
	#pragma omp parallel
	{
		seqan::AlignmentStats stats;

		#pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			computeAlignmentStats(stats, seqsh[i], seqsv[i], scoring_scheme);

			double alen_minus_gapopens =
				stats.alignmentLength - stats.numGapOpens;
			
			int len_seqh = seqan::length(seqan::source(seqsh[i]));
			int len_seqv = seqan::length(seqan::source(seqsv[i]));

			// only keep alignments that meet coverage and ani criteria
			if (std::max((alen_minus_gapopens / len_seqh),
						 (alen_minus_gapopens / len_seqv)) >= ratioScoreOverlap &&
				stats.alignmentIdentity >= debugThr)
			{
				dibella::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
				cks->score = (float)stats.alignmentIdentity / 100.0f;
				cks->passed = true;	// keep this
			}
		}
	}

	end_time = std::chrono::system_clock::now();
  	add_time("FA:ComputeStats + string_op",
			 (ms_t(end_time - start_time)).count());

	return;
}
