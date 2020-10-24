// Created by Saliya Ekanayake on 2019-09-03.

#include "../../include/pw/BandedAligner.hpp"

BandedAligner::BandedAligner(ScoringScheme scoring_scheme,
    int banded_half_width) :

    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    banded_half_width(banded_half_width){

	// int numThreads = 1;
	// #ifdef THREADED
	// #pragma omp parallel
	// {	
	// 	numThreads = omp_get_num_threads();
	// }
	// std::cout << "#threads " << numThreads << std::endl;
	// #endif

	// std::cout << "#threads (ctor ba) " << numThreads << std::endl;
	// tmp = new std::vector<int>[numThreads];
}

//template<typename TSequenceValue, typename TSpec>
//void SeedExtendXdrop<TSequenceValue, TSpec>::apply(
void BandedAligner::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
	ushort k,
    dibella::CommonKmers &cks, std::stringstream& ss) {

  auto start_pf_time = std::chrono::system_clock::now();

  seqan::Align<seqan::Dna5String> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), *seq_h);
  assignSource(row(align, 1), *seq_v);

  auto start_time = std::chrono::system_clock::now();
  int score = localAlignment(align, scoring_scheme, -banded_half_width, banded_half_width);
  auto end_time = std::chrono::system_clock::now();
  add_time("BA:local_alignment", (ms_t(end_time - start_time)).count());

  AlignmentInfo ai;
  start_time = std::chrono::system_clock::now();
  computeAlignmentStats(ai.stats, align, scoring_scheme);
  end_time = std::chrono::system_clock::now();
  add_time("BA:ComputeStats", (ms_t(end_time - start_time)).count());

  ai.seq_h_length = length(*seq_h);
  ai.seq_v_length = length(*seq_v);
  ai.seq_h_seed_length = (clippedEndPosition(row(align, 0)) - 1) - clippedBeginPosition(row(align, 0));
  ai.seq_v_seed_length = (clippedEndPosition(row(align, 1)) - 1) - clippedBeginPosition(row(align, 1));
  ai.seq_h_g_idx = g_col_idx;
  ai.seq_v_g_idx = g_row_idx;

//   /* Hard coding quality constraints for now */

//   // TODO - Saliya
//   // For now only keeps the largest alignment > 30% identity.
//   // Incorporate length coverage restrictions later.

// //  if (max_ai.stats.alignmentIdentity >= 30.0){
//   // alignments.push_back(ai);
// //  }

  double alen_minus_gapopens = (ai.stats.alignmentLength - ai.stats.numGapOpens) * 1.0;
  ss << g_col_idx << "," << g_row_idx << "," << ai.stats.alignmentIdentity
     << "," << ai.seq_h_length << "," << ai.seq_v_length
     << "," << ai.seq_h_seed_length  << "," << ai.seq_v_seed_length
     << "," << ai.stats.numGapOpens
     << "," << alen_minus_gapopens / ai.seq_h_length
     << "," << alen_minus_gapopens / ai.seq_v_length
	 << "," << cks.count
	 << std::endl;

  auto finish_pf_time = std::chrono::system_clock::now();
  add_time("BA:overall", (ms_t(finish_pf_time - start_pf_time)).count());
}

// Disabled for now (09/01/2020)
// @NOTE vectorized banded alignment sometimes hangs
void
BandedAligner::apply_batch
(
    seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
	seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
	uint64_t *lids,
	uint64_t col_offset,
	uint64_t row_offset,
    PSpMat<dibella::CommonKmers>::ref_tuples *mattuples,
    std::ofstream &lfs,
	ushort k,
    float ratioScoreOverlap,
	int debugThr
)
{
	// seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial> exec_policy;

	// int numThreads = 1;
	// #ifdef THREADED
	// #pragma omp parallel
    // {
    //   	numThreads = omp_get_num_threads();
    // }
	// #endif

	// uint64_t npairs = seqan::length(seqsh);
	// setNumThreads(exec_policy, numThreads);

	// lfs << "processing batch of size " << npairs << " with "
	// 	<< numThreads << " threads " << std::endl;

	// auto start_time = std::chrono::system_clock::now();
	
	// // alignment
	// localAlignment(exec_policy, seqsh, seqsv, scoring_scheme,
	// 			   -banded_half_width, banded_half_width);
	
	// auto end_time = std::chrono::system_clock::now();
  	// add_time("BA:local_alignment", (ms_t(end_time - start_time)).count());

	// start_time = std::chrono::system_clock::now();

	// lfs << "computing and writing stats" << std::endl;
	
	// // stats
	// #pragma omp parallel
	// {
	// 	seqan::AlignmentStats	stats;
	// 	std::stringstream		ss;

	// 	#pragma omp for schedule(static, 1000)
	// 	for (uint64_t i = 0; i < npairs; ++i)
	// 	{
	// 		computeAlignmentStats(stats, seqsh[i], seqsv[i], scoring_scheme);
			
	// 		double alen_minus_gapopens =
	// 			stats.alignmentLength - stats.numGapOpens;
	// 		int len_seqh = seqan::length(seqan::source(seqsh[i]));
	// 		int len_seqv = seqan::length(seqan::source(seqsv[i]));
	// 		ss << (col_offset + std::get<1>(mattuples[lids[i]])) << ","
	// 		   << (row_offset + std::get<0>(mattuples[lids[i]]))  << ","
	// 		   << stats.alignmentIdentity << ","
	// 		   << len_seqh << ","
	// 		   << len_seqv << ","
	// 		   << (clippedEndPosition(seqsh[i]) -
	// 			   clippedBeginPosition(seqsh[i]) - 1) << ","
	// 		   << (clippedEndPosition(seqsv[i]) -
	// 			   clippedBeginPosition(seqsv[i]) - 1) << ","
	// 		   << stats.numGapOpens << ","
	// 		   << alen_minus_gapopens / len_seqh << ","
	// 		   << alen_minus_gapopens / len_seqv
	// 		   << "\n";
	// 	}

	// 	#pragma omp critical
	// 	{
	// 		afs << ss.str();
	// 		afs.flush();
	// 	}
	// }

	// end_time = std::chrono::system_clock::now();
  	// add_time("BA:ComputeStats + string_op",
	// 		 (ms_t(end_time - start_time)).count());

	// return;
}
