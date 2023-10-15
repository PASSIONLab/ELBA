// Created by Saliya Ekanayake on 2019-07-05 and modified by Giulia Guidi on 09/01/2020.

#include "../../include/pw/SeedExtendXdrop.hpp"
#include <unordered_set>

uint minOverlapLen = 5000;

void SeedExtendXdrop::PostAlignDecision(const AlignmentInfo& ai, bool& passed, float& ratioScoreOverlap,
	int& dir, int& dirT, int& sfx, int& sfxT, uint32_t& overlap, const bool noAlign, std::vector<int64_t>& ContainedSeqMyThread)
{
	auto maxseed = ai.seed;	// returns a seqan:Seed object

	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the query (vertical/horizonral direction)
	// these four return seqan:Tposition objects
	int begpV = beginPositionV(maxseed);
	int endpV = endPositionV  (maxseed);
	int begpH = beginPositionH(maxseed);
	int endpH = endPositionH  (maxseed);

	unsigned short int overlapLenH = ai.seq_h_seed_length;
	unsigned short int overlapLenV = ai.seq_v_seed_length;

	unsigned short int rlenH = ai.seq_h_length;
	unsigned short int rlenV = ai.seq_v_length;

	unsigned short int minLeft  = min(begpV, begpH);
	unsigned short int minRight = min(rlenV - endpV, rlenH - endpH);

	int64_t seqV = ai.seq_v_g_idx;
	int64_t seqH = ai.seq_h_g_idx;

	overlap = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

	if((seqV == 5 && seqH == 99) || 
		(seqV == 5 && seqH == 100) || 
			(seqV == 5 && seqH == 184) ||
				(seqV == 24 && seqH == 40))
	{
		std::cout << seqV+1 << "\t" << seqH+1 << "\t" << ai.rc << "\t" << begpV 
		<< "\t" << endpV << "\t" << begpH << "\t" << endpH  << "\t" << ai.xscore << "\t" << overlap << std::endl;
	}

#ifndef FIXEDTHR
	float myThr = (1 - DELTACHERNOFF) * (ratioScoreOverlap * (float)overlap);

	// Contained overlaps removed for now, reintroduce them later
	// @GGGG-TODO: identify chimeric sequences
	bool contained = false;
	bool chimeric  = false;

	if (begpV <= begpH && (rlenV - endpV) <= (rlenH - endpH))
	{
	    ContainedSeqMyThread.push_back(seqV);
	    contained = true;
	}
	else if (begpV >= begpH && (rlenV - endpV) >= (rlenH - endpH))
	{
	    ContainedSeqMyThread.push_back(seqH);
	    contained = true;
	}
	else if (!noAlign)
	{
	    passed = ((float)ai.xscore >= myThr && overlap >= minOverlapLen);

	    if (passed)
	    {
	        if (begpV > begpH)
	        {
	            dir  = ai.rc? 0 : 1;
	            dirT = ai.rc? 0 : 2;
	            sfx  = ((rlenH - endpH) - (rlenV - endpV));
	            sfxT = begpV - begpH;
	        }
	        else
	        {
	            dir  = ai.rc? 3 : 2;
	            dirT = ai.rc? 3 : 1;
	            sfx  = begpH - begpV;
	            sfxT = ((rlenV - endpV) - (rlenH - endpH));
	        }
	    }
	}

#else
	if(ai.xscore >= FIXEDTHR)
		passed = true;
#endif
}

SeedExtendXdrop::SeedExtendXdrop(
    ScoringScheme scoring_scheme,
    ushort seed_length, int xdrop, int seed_count):
    PairwiseFunction(),
    scoring_scheme(scoring_scheme),
    seed_length(seed_length), xdrop(xdrop), seed_count(seed_count){
}

void SeedExtendXdrop::apply(
    uint64_t l_col_idx, uint64_t g_col_idx,
    uint64_t l_row_idx, uint64_t g_row_idx,
    seqan::Dna5String *seqH, seqan::Dna5String *seqV, ushort k,
    elba::CommonKmers &cks, std::stringstream& ss)
{
    // ...
}

// @NOTE This is hard-coded to the number of seeds being <= 2
void
SeedExtendXdrop::apply_batch
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
	std::vector<int64_t>& ContainedSeqPerBatch,
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

	lfs << "processing batch of size " << npairs << " with " << numThreads << " threads " << std::endl;

	// for multiple seeds we store the seed with the highest identity
	AlignmentInfo *ai = new AlignmentInfo[npairs];
	std::pair<ushort, ushort> *seedlens = new std::pair<ushort, ushort>[npairs];

	bool *strands = new bool[npairs];
	int  *xscores = new int[npairs];
	TSeed  *seeds = new TSeed[npairs];
		
//	int runseqan = 0;


	/* GGGG: seed_count is hardcoded here (2) */
	for(int count = 0; count < seed_count; ++count)
	{
		auto start_time = std::chrono::system_clock::now();

		seqan::StringSet<seqan::Dna5String> seqsh_ex;
		seqan::StringSet<seqan::Dna5String> seqsv_ex;
		resize(seqsh_ex, npairs, seqan::Exact{});
		resize(seqsv_ex, npairs, seqan::Exact{});

		// extend the current seed and form a new gaps object
		#pragma omp parallel for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

			// In KmerIntersectSR.hpp we have (where res == cks):
			// 	res.first.first 	// Kmer 1 on argA
			// 	res.first.second 	// Kmer 1 on argB
			// 	res.second.first 	// Kmer 2 on argA
			// 	res.second.second 	// Kmer 2 on argB

			// argA (see KmerIntersectSR.hpp) == row == seqV
			ushort LocalSeedVOffset = (count == 0) ? cks->first.first : cks->second.first;
			// argB (see KmerIntersectSR.hpp) == col == seqH
			ushort LocalSeedHOffset = (count == 0) ? cks->first.second : cks->second.second;

			seqan::Dna5String seedV; // seed on arg1 == row == seqV
			seqan::Dna5String seedH; // seed on arg2 == col == seqH

			auto start_time = std::chrono::system_clock::now();
			auto end_time   = std::chrono::system_clock::now();

			// seed creation params are:
			// horizontal seed start offset, vertical seed start offset, length
			TSeed seed(LocalSeedHOffset, LocalSeedVOffset, seed_length);

			seedV = infix(seqsv[i], beginPositionV(seed), endPositionV(seed)); // seed on argA == row == seqV
			seedH = infix(seqsh[i], beginPositionH(seed), endPositionH(seed)); // seed on argB == col == seqH

			seqan::Dna5StringReverseComplement twin(seedH);

			if(twin == seedV)
			{
				strands[i] = true;
				seqan::Dna5String twinseqH = seqsh[i];
				seqan::Dna5StringReverseComplement twinRead(twinseqH);
				LocalSeedHOffset = length(twinseqH) - LocalSeedHOffset - seed_length;

				setBeginPositionH(seed, LocalSeedHOffset);
				setBeginPositionV(seed, LocalSeedVOffset);
				setEndPositionH(seed, LocalSeedHOffset + seed_length);
				setEndPositionV(seed, LocalSeedVOffset + seed_length);

				if(!noAlign)
				{
					/* Perform match extension */
					start_time = std::chrono::system_clock::now();
	//				runseqan++;
					xscores[i] = extendSeed(seed, twinRead, seqsv[i], seqan::EXTEND_BOTH, scoring_scheme,
							xdrop, (int)k,
							seqan::GappedXDrop());

					end_time = std::chrono::system_clock::now();
					add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());
				}
				else
				{
					xscores[i] = 0;
				}
			}
			else if(seedH == seedV)
			{

				strands[i] = false;
				if(!noAlign)
				{
					start_time = std::chrono::system_clock::now();
	//				runseqan++;
					xscores[i] = extendSeed(seed, seqsh[i], seqsv[i], seqan::EXTEND_BOTH, scoring_scheme,
							xdrop, (int)k,
							seqan::GappedXDrop());
					end_time = std::chrono::system_clock::now();
					add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

				}
				else
				{
					xscores[i] = 0;
				}
			}

			seeds[i] = seed;
			seedlens[i].first  = static_cast<ushort>(seed._endPositionH -
													 seed._beginPositionH);
			seedlens[i].second = static_cast<ushort>(seed._endPositionV -
													 seed._beginPositionV);
		}

		auto end_time = std::chrono::system_clock::now();
    		add_time("XA:ExtendSeed", (ms_t(end_time - start_time)).count());

		start_time = std::chrono::system_clock::now();

//		std::cout << "runseqan " << runseqan << std::endl;

		// Compute stats
		if (count == 0)	// GGGG: save the results from the first alignment (k-mer 1)
		{
		#pragma omp parallel for
			for (uint64_t i = 0; i < npairs; ++i)
			{
      
				ai[i].xscore = xscores[i];
				ai[i].rc     = strands[i];
				ai[i].seed   =   seeds[i];

				ai[i].seq_h_length = seqan::length(seqsh[i]);
				ai[i].seq_v_length = seqan::length(seqsv[i]);

				ai[i].seq_h_seed_length = seedlens[i].first;
				ai[i].seq_v_seed_length = seedlens[i].second;

				// GGGG: global idx over here to use in the FullDistVect for removing contained vertices/seqs
				ai[i].seq_h_g_idx = col_offset + std::get<1>(mattuples[lids[i]]);
    			ai[i].seq_v_g_idx = row_offset + std::get<0>(mattuples[lids[i]]);
			}
		}
		else // GGGG: if the results from the second seed if better, then substitute it in ai
		{
		#pragma omp parallel for
	 	for (uint64_t i = 0; i < npairs; ++i)
		 	{
				if (xscores[i] > ai[i].xscore) // GGGG: if new seed is better, then save the new one
		 		//  GGGG: but are we sure the score standalons is the right metric?
		 		{
		 			ai[i].xscore = xscores[i];
		 			ai[i].rc     = strands[i];
		 			ai[i].seed   =   seeds[i];
		 			ai[i].seq_h_seed_length = seedlens[i].first;
	 			ai[i].seq_v_seed_length = seedlens[i].second;
		 		}
		 	}
		 }

		end_time = std::chrono::system_clock::now();
    	add_time("XA:ComputeStats", (ms_t(end_time - start_time)).count());
	}

	delete [] seedlens;
	delete [] xscores;
	delete [] strands;
	delete [] seeds;

	auto start_time = std::chrono::system_clock::now();

	std::vector<std::vector<int64_t>> ContainedSeqPerThread(numThreads);

	// TODO@GGGG: reproduce and fix segfault
	// for(int t = 0; t < numThreads; t++)
	// 	ContainedSeqPerThread[t].resize(std::ceil(nreads/numThreads));

	// Dump alignment info
	#pragma omp parallel
	{
	  #pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			// Only keep alignments that meet BELLA criteria
			bool passed = false;
			int tid = omp_get_thread_num();

			elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

			// GGGG: ai stores global idx to to store in ContainedSeqPerBatch
			// GGGG: in PostAlignDecision() we can mark as contained sequences as removable in ContainedSeqPerBatch and their local contained edges
			// GGGG: ContainedSeqPerBatch global indexes of contained sequences

			PostAlignDecision(ai[i], passed, ratioScoreOverlap, cks->dir, cks->dirT, cks->sfx, cks->sfxT, cks->overlap, noAlign, ContainedSeqPerThread[tid]);

			if (passed) // GGGG: this should probably decide which seed, rather than earlier?
			{
				// GGGG: store updated seed start/end position in the CommonKmers pairs (the semantics of these pairs change wrt the original semantics but that's okay)
				cks->first.first   = beginPositionV(ai[i].seed); 	// start on vertical sequence
				cks->first.second  = endPositionV(ai[i].seed); 		// end on vertical sequence
				cks->second.first  = beginPositionH(ai[i].seed);	// start on horizonal sequence
				cks->second.second = endPositionH(ai[i].seed);		// end on horizonal sequence

				cks->lenv 	= ai[i].seq_v_length;
				cks->lenh 	= ai[i].seq_h_length;
				cks->score  = ai[i].xscore;
                cks->rc     = ai[i].rc;
				cks->passed = passed;	// keep this
			}
		}
	}

	int readcount = 0;
	for(int t = 0; t < numThreads; ++t)
	{
		readcount += ContainedSeqPerThread[t].size();
	}

	unsigned int readssofar = 0;
	ContainedSeqPerBatch.resize(readcount);

	// Concatenate per-thread result
	for(int t = 0; t < numThreads; ++t)
	{
		copy(ContainedSeqPerThread[t].begin(), ContainedSeqPerThread[t].end(), ContainedSeqPerBatch.begin() + readssofar);
		readssofar += ContainedSeqPerThread[t].size();
	}

	auto end_time = std::chrono::system_clock::now();
  	add_time("XA:StringOp", (ms_t(end_time - start_time)).count());

	delete [] ai;

	return;
}
