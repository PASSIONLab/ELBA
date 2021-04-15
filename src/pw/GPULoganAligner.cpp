/* Created by Saliya Ekanayake on 2019-07-05 and modified by Giulia Guidi on 4/14/2021. */

// @TODO: import load balancer from BELLA

#include "../../include/pw/GPULoganAligner.hpp"

// uint minOverlapLen = 10000;

// void GPULoganAligner::PostAlignDecision(const AlignmentInfo& ai, bool& passed, float& ratioScoreOverlap, 
// 	uint32_t& overhang, uint32_t& overhangT, uint32_t& overlap, const bool noAlign)
// {
// 	auto maxseed = ai.seed;	// returns a seqan:Seed object

// 	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the query (vertical/horizonral direction)
// 	// these four return seqan:Tposition objects
// 	int begpV = beginPositionV(maxseed);
// 	int endpV = endPositionV  (maxseed);
// 	int begpH = beginPositionH(maxseed);
// 	int endpH = endPositionH  (maxseed);

// 	unsigned short int overlapLenH = ai.seq_h_seed_length;
// 	unsigned short int overlapLenV = ai.seq_v_seed_length;

// 	unsigned short int rlenH = ai.seq_h_length;
// 	unsigned short int rlenV = ai.seq_v_length;

// 	unsigned short int minLeft  = min(begpV, begpH);
// 	unsigned short int minRight = min(rlenV - endpV, rlenH - endpH);

// 	overlap = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

// #ifndef FIXEDTHR
// 	float myThr = (1 - DELTACHERNOFF) * (ratioScoreOverlap * (float)overlap);

// 	// Contained overlaps removed for now, reintroduce them later
// 	bool contained = false;
// 	bool chimeric  = false; 

// 	// seqH is column entry and seqV is row entry, for each column, we iterate over seqHs
// 	int seqH = ai.seq_v_g_idx, seqV = ai.seq_h_g_idx;
	
// 	// Reserve length/position if rc [x]
// 	if(ai.rc)
// 	{
// 		uint tmp = begpV;
// 		begpV = rlenV - endpV;
// 		endpV = rlenV - tmp;
// 	}

// 	if((begpH == 0 & rlenH-endpH == 0) || (begpV == 0 & rlenV-endpV == 0))
// 		contained = true;
	
// 	if(!contained)
// 	{
// 		// If noAlign is false, set passed to false if the score isn't good enough
// 		if(!noAlign)
// 		{
// 			if((float)ai.xscore < myThr || overlap < minOverlapLen) passed = false;
// 			else passed = true;
// 		}

// 		if(passed)
// 		{
// 			uint32_t direction, directionT;
// 			uint32_t suffix, suffixT;

// 			// !reverse complement
// 			if(!ai.rc)
// 			{
// 				if(begpH > begpV)
// 				{
// 					direction  = 1;
// 					directionT = 2;

// 					suffix  = rlenV - endpV;
// 					suffixT = begpH;
// 				}	
// 				else
// 				{
// 					direction  = 2;
// 					directionT = 1;

// 					suffix  = rlenH - endpH;
// 					suffixT = begpV;
// 				} 
// 			}
// 			else
// 			{
// 				if((begpV > 0) & (begpH > 0) & (rlenV-endpV == 0) & (rlenV-endpV == 0))
// 				{
// 					direction  = 0;
// 					directionT = 0;

// 					suffix  = begpV; // seqV == 2
// 					suffixT = begpH; // seqV == 2			
// 				}
// 				else
// 				{
// 					direction  = 3;
// 					directionT = 3;

// 					suffix  = rlenV - endpV; // seqV == 1, seqH == 2	
// 					suffixT = rlenH - endpH; // seqV == 1, seqH == 2		
// 				}
// 			}
// 			overhang  = suffix  << 2 | direction;
// 			overhangT = suffixT << 2 | directionT;
// 		} // if(passed)
// 	} // if(!contained)
		
// #else
// 	if(ai.xscore >= FIXEDTHR)
// 		passed = true;
// #endif
// }

GPULoganAligner::GPULoganAligner(
    ScoringScheme scoring_scheme,
    ushort seed_length, int xdrop, int seed_count):
    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    seed_length(seed_length), xdrop(xdrop), seed_count(seed_count){
}

void GPULoganAligner::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seqH, seqan::Dna5String *seqV, ushort k,
    dibella::CommonKmers &cks, std::stringstream& ss)
{
    // ...
}

// @NOTE This is hard-coded to the number of seeds being <= 2
void
GPULoganAligner::apply_batch
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
    float ratioScoreOverlap, // GGGG: this is my ratioScoreOverlap variable change name later
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

	// for multiple seeds we store the seed with the highest identity
	AlignmentInfo *ai = new AlignmentInfo[npairs];
	std::pair<ushort, ushort> *seedlens = new std::pair<ushort, ushort>[npairs];

	// bool *strands = new bool[npairs]; // Directly in loganResult
	// int  *xscores = new int[npairs];
	// TSeed  *seeds = new TSeed[npairs];

	std::vector<string> seqHs;
	std::vector<string> seqVs;
	std::vector<SeedL>  seeds;
	std::vector<loganResult> xscores;

	/* GGGG: seed_count is hardcoded here (2) */
	for(int count = 0; count < seed_count; ++count)
	{
		auto start_time = std::chrono::system_clock::now();

		// seqan::StringSet<seqan::Gaps<seqan::Dna5String>> seqsh_ex;
		// seqan::StringSet<seqan::Gaps<seqan::Dna5String>> seqsv_ex;
		// resize(seqsh_ex, npairs, seqan::Exact{});
		// resize(seqsv_ex, npairs, seqan::Exact{});
	
		// @GGGG: keep the order for the post alignment evaluation (this might not be feasible in this case? measure slowdown)
		// #pragma omp parallel for 
		for(IT j = start; j < end; ++j) // I acculate sequences for GPU batch alignment
			for (uint64_t i = 0; i < npairs; ++i)
			{
				dibella::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

			#ifdef TWOSEED
				ushort LocalSeedVOffset =
					(count == 0) ? cks->first.first : cks->second.first;
				ushort LocalSeedHOffset =
					(count == 0) ? cks->first.second : cks->second.second;
			#else
				// GGGG: TODO check reverse complement
				ushort LocalSeedVOffset = cks.pos[0].first;
				ushort LocalSeedHOffset = cks.pos[0].second;
			#endif

				seqan::Dna5String seedH;
				seqan::Dna5String seedV;

				auto start_time = std::chrono::system_clock::now();
				auto end_time   = std::chrono::system_clock::now();

				// @GGGG-TODO: remove this stuff and only use string
				seedH = infix(seqan::source(seqsh[i]), beginPositionH(seed), endPositionH(seed));
				seedV = infix(seqan::source(seqsv[i]), beginPositionV(seed), endPositionV(seed)); 

				seqan::Dna5StringReverseComplement twin(seedH);

			#ifdef STATS
				seqan::Align<seqan::Dna5String> align;
				resize(rows(align), 2);
			#endif

				std::string seqV = seqsv[i];
				std::string seqH = seqsh[i];

				logan::loganResult localRes; // not sure we need logan::

				if(twin == seedV)
				{
					localRes.strand = true;

					std::string twinseqH(seqH);

					// @GGGG-TODO: remove this stuff and only use string
					seqan::Dna5String twinseqH = seqan::source(seqsh[i]);
					seqan::Dna5StringReverseComplement twinRead(twinseqH);
					LocalSeedHOffset = length(twinseqH) - LocalSeedHOffset - seed_length;

					logan::SeedL seed(LocalSeedHOffset, LocalSeedVOffset, LocalSeedHOffset + seed_length, LocalSeedVOffset + seed_length);

					if(!noAlign) 
					{
						// GGGG: here only accumulate stuff for the GPUs, don't perform alignment
						seeds.push_back(seed);
						seqVs.push_back(seqV);
						seqHs.push_back(twinseqH);

						xscores.push_back(localRes);

					#ifdef STATS
						assignSource(row(align, 0), infix(twinRead, beginPositionH(seed),
														endPositionH(seed)));
						assignSource(row(align, 1), infix(*seqV, beginPositionV(seed),
														endPositionV(seed)));
					#endif
					}
					else // @GGGG-TODO: the logic is a bit different than BELLA -- double check and clean unused code
					{
						localRes.xscore = 0;

						seeds.push_back(seed);
						seqVs.push_back(seqV);
						seqHs.push_back(twinseqH);

						xscores.push_back(localRes);
					}
				}
				else
				{
					localRes.strand = false;

					// @GGGG-TODO: remove this stuff and only use string
					seqan::Dna5String twinseqH = seqan::source(seqsh[i]);
					seqan::Dna5StringReverseComplement twinRead(twinseqH);
					LocalSeedHOffset = length(twinseqH) - LocalSeedHOffset - seed_length;

					logan::SeedL seed(LocalSeedHOffset, LocalSeedVOffset, LocalSeedHOffset + seed_length, LocalSeedVOffset + seed_length);

					if(!noAlign) 
					{
						// GGGG: here only accumulate stuff for the GPUs, don't perform alignment
						seeds.push_back(seed);
						seqVs.push_back(seqV);
						seqHs.push_back(twinseqH);

						xscores.push_back(localRes);
						
					#ifdef STATS
						assignSource(row(align, 0), infix(twinRead, beginPositionH(seed),
														endPositionH(seed)));
						assignSource(row(align, 1), infix(*seqV, beginPositionV(seed),
														endPositionV(seed)));
					#endif
					}
					else // @GGGG-TODO: the logic is a bit different than BELLA -- double check and clean unused code
					{
						localRes.xscore = 0;

						seeds.push_back(seed);
						seqVs.push_back(seqV);
						seqHs.push_back(twinseqH);

						xscores.push_back(localRes);
					}
				}

			#ifdef STATS
				xscores[i] = extendSeed(seed, seqan::source(seqsh[i]), seqan::source(seqsv[i]),
						seqan::EXTEND_BOTH, scoring_scheme,
						xdrop, (int)k, seqan::GappedXDrop());
				assignSource(seqsh_ex[i],
							infix(seqan::source(seqsh[i]),
								beginPositionH(seed), endPositionH(seed)));
				assignSource(seqsv_ex[i],
							infix(seqan::source(seqsv[i]),
								beginPositionV(seed), endPositionV(seed)));
			#endif
				seeds[i] = seed;
				seedlens[i].first  = static_cast<ushort>(seed._endPositionH -
														seed._beginPositionH);
				seedlens[i].second = static_cast<ushort>(seed._endPositionV -
														seed._beginPositionV);
			}

		auto end_time = std::chrono::system_clock::now();
    	add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

	#ifdef STATS
		start_time = std::chrono::system_clock::now();
		// alignment
		globalAlignment(exec_policy, seqsh_ex, seqsv_ex, scoring_scheme);
		
		end_time = std::chrono::system_clock::now();
    	add_time("XA:global_alignment", (ms_t(end_time - start_time)).count());
	#endif
		start_time = std::chrono::system_clock::now();
		
		// Compute stats
		if (count == 0)	// overwrite in the first seed
		{
		#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
			#ifdef STATS
				computeAlignmentStats(ai[i].stats, seqsh_ex[i], seqsv_ex[i],
									  scoring_scheme);
			#endif
				ai[i].xscore = xscores[i];
				ai[i].rc     = strands[i];
				ai[i].seed   =   seeds[i];

				ai[i].seq_h_length = seqan::length(seqan::source(seqsh[i]));
				ai[i].seq_v_length = seqan::length(seqan::source(seqsv[i]));

				ai[i].seq_h_seed_length = seedlens[i].first;
				ai[i].seq_v_seed_length = seedlens[i].second;

				ai[i].seq_h_g_idx = col_offset + std::get<1>(mattuples[lids[i]]);
    			ai[i].seq_v_g_idx = row_offset + std::get<0>(mattuples[lids[i]]);
			}
		}
		else
		{
		#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{

			#ifdef STATS
				seqan::AlignmentStats stats;
				computeAlignmentStats(stats, seqsh_ex[i], seqsv_ex[i],
									  scoring_scheme);
		
				if (stats.alignmentIdentity > ai[i].stats.alignmentIdentity)
				{
					ai[i].stats				= stats;
					ai[i].seq_h_seed_length = seedlens[i].first;
					ai[i].seq_v_seed_length = seedlens[i].second;
				}
			#else
				if (xscores[i] > ai[i].xscore) // GGGG: TODO double check this logic with fresh neurons
				{
					ai[i].xscore = xscores[i];
					ai[i].rc     = strands[i];
					ai[i].seed   =   seeds[i];
					ai[i].seq_h_seed_length = seedlens[i].first;
					ai[i].seq_v_seed_length = seedlens[i].second;
				}
			#endif
			}
		}

		end_time = std::chrono::system_clock::now();
    	add_time("XA:ComputeStats", (ms_t(end_time - start_time)).count());
	}


	delete [] seedlens;
	delete [] xscores;
	delete [] strands;

	auto start_time = std::chrono::system_clock::now();
	
	// Dump alignment info
	#pragma omp parallel
	{
	    #pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			// Only keep alignments that meet BELLA criteria
			bool passed = false;

			dibella::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
			PostAlignDecision(ai[i], passed, ratioScoreOverlap, cks->overhang, cks->overhangT, cks->overlap, noAlign);

			if (passed)
			{
				// GGGG: store updated seed start/end position in the CommonKmers pairs (the semantics of these pairs change wrt the original semantics but that's okay)
				cks->first.first   = beginPositionV(ai[i].seed); 	// start on vertical sequence
				cks->first.second  = endPositionV(ai[i].seed); 		// end on vertical sequence
				cks->second.first  = beginPositionH(ai[i].seed);	// start on horizonal sequence
				cks->second.second = endPositionH(ai[i].seed);		// end on horizonal sequence

				cks->lenv 	= ai[i].seq_v_length;
				cks->lenh 	= ai[i].seq_h_length;
				cks->score  = ai[i].xscore;
				cks->passed = passed;	// keep this
			}
		}
	}

	auto end_time = std::chrono::system_clock::now();
  	add_time("XA:StringOp",
			 (ms_t(end_time - start_time)).count());

	delete [] ai;

	return;
}
