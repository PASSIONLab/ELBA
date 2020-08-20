// Created by Saliya Ekanayake on 2019-09-03.

#include "../../include/pw/FullAligner.hpp"

FullAligner::FullAligner(ScoringScheme scoring_scheme) :
    PairwiseFunction(),
    scoring_scheme(scoring_scheme)
	{}

//template<typename TSequenceValue, typename TSpec>
//void SeedExtendXdrop<TSequenceValue, TSpec>::apply(
void FullAligner::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
    dibella::CommonKmers &cks, std::stringstream &ss) {

  seqan::Align<seqan::Dna5String> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), *seq_h);
  assignSource(row(align, 1), *seq_v);

  auto start_time = std::chrono::system_clock::now();
  int score = localAlignment(align, scoring_scheme, seqan::DynamicGaps());
  auto end_time = std::chrono::system_clock::now();
  add_time("FA:local_alignment", (ms_t(end_time - start_time)).count());

  AlignmentInfo ai;
  start_time = std::chrono::system_clock::now();
  computeAlignmentStats(ai.stats, align, scoring_scheme);
  end_time = std::chrono::system_clock::now();
  add_time("FA:compute_stats", (ms_t(end_time - start_time)).count());

  ai.seq_h_length = length(*seq_h);
  ai.seq_v_length = length(*seq_v);
  ai.seq_h_seed_length = (clippedEndPosition(row(align, 0)) - 1) -
                         clippedBeginPosition(row(align, 0));
  ai.seq_v_seed_length = (clippedEndPosition(row(align, 1)) - 1) -
                         clippedBeginPosition(row(align, 1));
  ai.seq_h_g_idx = g_col_idx;
  ai.seq_v_g_idx = g_row_idx;

  /* Hard coding quality constraints for now */

  // TODO - Saliya
  // For now only keeps the largest alignment > 30% identity.
  // Incorporate length coverage restrictions later.

  // TODO - Keep a counter
//  if (max_ai.stats.alignmentIdentity >= 30.0){
//  alignments.push_back(ai);
//  }

  start_time = std::chrono::system_clock::now();
  double alen_minus_gapopens = (ai.stats.alignmentLength - ai.stats.numGapOpens) * 1.0;
  ss << g_col_idx << "," << g_row_idx << "," << ai.stats.alignmentIdentity
     << "," << ai.seq_h_length << "," << ai.seq_v_length
     << "," << ai.seq_h_seed_length << "," << ai.seq_v_seed_length
     << "," << ai.stats.numGapOpens
     << "," << alen_minus_gapopens / ai.seq_h_length
     << "," << alen_minus_gapopens / ai.seq_v_length
	 << "," << cks.count
	 << std::endl;
  end_time = std::chrono::system_clock::now();
  add_time("FA:string_op", (ms_t(end_time - start_time)).count());
}



void
FullAligner::apply_batch
(
    seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
	seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
	uint64_t *lids,
	uint64_t col_offset,
	uint64_t row_offset,
	PSpMat<dibella::CommonKmers>::Tuples &mattuples,
	std::ofstream &afs,
	std::ofstream &lfs
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
		seqan::AlignmentStats	stats;
		std::stringstream		ss;

		#pragma omp for schedule(static, 1000)
		for (uint64_t i = 0; i < npairs; ++i)
		{
			computeAlignmentStats(stats, seqsh[i], seqsv[i], scoring_scheme);
			
			double alen_minus_gapopens =
				stats.alignmentLength - stats.numGapOpens;
			int len_seqh = seqan::length(seqan::source(seqsh[i]));
			int len_seqv = seqan::length(seqan::source(seqsv[i]));
			ss << (col_offset + mattuples.colindex(lids[i])) << ","
			   << (row_offset + mattuples.rowindex(lids[i]))  << ","
			   << stats.alignmentIdentity << ","
			   << len_seqh << ","
			   << len_seqv << ","
			   << (clippedEndPosition(seqsh[i]) -
				   clippedBeginPosition(seqsh[i]) - 1) << ","
			   << (clippedEndPosition(seqsv[i]) -
				   clippedBeginPosition(seqsv[i]) - 1) << ","
			   << stats.numGapOpens << ","
			   << alen_minus_gapopens / len_seqh << ","
			   << alen_minus_gapopens / len_seqv
			   << "\n";
		}

		#pragma omp critical
		{
			afs << ss.str();
			afs.flush();
		}
	}

	end_time = std::chrono::system_clock::now();
  	add_time("FA:compute_stats + string_op",
			 (ms_t(end_time - start_time)).count());

	return;
}
