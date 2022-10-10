/* Created by Saliya Ekanayake on 2019-07-05 and modified by Giulia Guidi on 4/14/2021. */

#include "../../include/pw/GPULoganAligner.hpp"
#include "../../LoganGPU/RunLoganAligner.hpp" 	// Call to aligner

uint minOverlapLenL = 5000;
uint maxOverhangL = 25;

char 
complementbase(char n)
{   
    switch(n)
    {   
    case 'A':
        return 'T';
    case 'T':
        return 'A';
    case 'G':
        return 'C';
    case 'C':
        return 'G';
    }   
    assert(false);
    return ' ';
} 

std::string
reversecomplement(const std::string& seq) {

	std::string cpyseq = seq;
	std::reverse(cpyseq.begin(), cpyseq.end());

	std::transform(
		std::begin(cpyseq),
		std::end  (cpyseq),
		std::begin(cpyseq),
	complementbase);

	return cpyseq;
}

void GPULoganAligner::PostAlignDecision(const LoganAlignmentInfo& ai, bool& passed, float& ratioScoreOveralap,
        int& dir, int& dirT, int& sfx, int& sfxT, uint32_t& overlap, const bool noAlign, std::vector<int64_t>& ContainedSeqPerProc);
{
	// {begin/end}Position{V/H}: Returns the begin/end position of the seed in the seqVs (vertical/horizonral direction)
	// these four return seqan:Tposition objects
	int begpV = ai.begSeedV;
	int endpV = ai.endSeedV;
	int begpH = ai.begSeedH;
	int endpH = ai.endSeedH;

	unsigned short int overlapLenH = ai.seq_h_seed_length;
	unsigned short int overlapLenV = ai.seq_v_seed_length;

	unsigned short int rlenH = ai.seq_h_length;
	unsigned short int rlenV = ai.seq_v_length;

	unsigned short int minLeft  = min(begpV, begpH);
	unsigned short int minRight = min(rlenV - endpV, rlenH - endpH);

	int64_t seqV = ai.seq_v_g_idx;
	int64_t seqH = ai.seq_h_g_idx;

	overlap = minLeft + minRight + (overlapLenV + overlapLenH) / 2;

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

	// // Reserve length/position if rc [x]
	// if(ai.rc)
	// {
	// 	uint tmp = begpH;
	// 	begpH = rlenH - endpH;
	// 	endpH = rlenH - tmp;
	// }

    // if (begpV > begpH && (rlenV - endpV) > (rlenH - endpH))
	// {
	// 	ContainedSeqMyThread.push_back(seqH); // Push back global index
	// 	contained = true;
	// }
    // else if (begpH > begpV && (rlenH - endpH) > (rlenV - endpV))
	// {
	// 	ContainedSeqMyThread.push_back(seqV); // Push back global index
	// 	contained = true;
	// }

	// if(!contained)
	// {
	// 	// If noAlign is false, set passed to false if the score isn't good enough
	// 	if(!noAlign)
	// 	{
	// 		if((float)ai.xscore < myThr || overlap < minOverlapLenL) passed = false;
	// 		else passed = true;
	// 	}

	// 	if(passed)
	// 	{
	// 		uint32_t direction, directionT;
	// 		uint32_t suffix, suffixT;

	// 		// !reverse complement
	// 		if(!ai.rc)
	// 		{
	// 			if(begpV > begpH)
	// 			{
	// 				direction  = 1;
	// 				directionT = 2;

	// 				suffix = rlenH - endpH;
	// 				suffixT = begpV;
	// 			}
	// 			else
	// 			{
	// 				direction  = 2;
	// 				directionT = 1;

	// 				suffix = begpH;
	// 				suffixT = rlenV - endpV;
	// 			}
	// 		}
	// 		else
	// 		{
	// 			if((begpV > 0) && (begpH > 0) && (rlenV-endpV == 0) && (rlenH-endpH == 0))
	// 			{
	// 				direction  = 0;
	// 				directionT = 0;

	// 				suffix  = begpH;
	// 				suffixT = begpV;
	// 			}
	// 			else
	// 			{
	// 				direction  = 3;
	// 				directionT = 3;

	// 				suffix  = rlenH - endpH;
	// 				suffixT = rlenV - endpV;
	// 			}
	// 		}
	// 		overhang  = suffix  << 2 | direction;
	// 		overhangT = suffixT << 2 | directionT;
	// 	} // if(passed)
	// } // if(!contained)

#else
	if(ai.xscore >= FIXEDTHR)
		passed = true;
#endif
}

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
    elba::CommonKmers &cks, std::stringstream& ss)
{
    // ...
}

// @NOTE This is hard-coded to the number of seeds being <= 2
void
GPULoganAligner::apply_batch
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
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	uint64_t npairs = seqan::length(seqsh);
	
	lfs << "processing batch of size " << npairs << " with " << numThreads << " threads " << std::endl;

	// for multiple seeds we store the seed with the highest identity
	LoganAlignmentInfo *ai = new LoganAlignmentInfo[npairs];
	std::pair<ushort, ushort> *seedlens = new std::pair<ushort, ushort>[npairs];

	std::vector<string> seqHs;
	std::vector<string> seqVs;

	std::vector<SeedInterface> seeds;
	std::vector<LoganResult> xscores;

	// bool *strands = new bool[npairs];
	// int  *xscores = new int[npairs];
	// TSeed  *seeds = new TSeed[npairs];

	/* GGGG: seed_count is hardcoded here (2) */
	for(int count = 0; count < seed_count; ++count)
	{
		auto start_time = std::chrono::system_clock::now();
	
		// @GGGG: keep the order for the post alignment evaluation (measure slowdown)
		// #pragma omp parallel for 
		for (uint64_t i = 0; i < npairs; ++i) // I acculate sequences for GPU batch alignment
		{
			// init result
			LoganResult localRes; 

			// Get seed location
			elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);

			// In KmerIntersectSR.hpp we have (where res == cks):
			// 	res.first.first 	= arg1.first.first;		// Kmer 1 on argA
			// 	res.first.second 	= arg1.first.second;	// Kmer 2 on argA
			// 	res.second.first 	= arg2.first.first;		// Kmer 1 on argB
			// 	res.second.second 	= arg2.first.second;	// Kmer 2 on argB

			// argA (see KmerIntersectSR.hpp) == row == seqV
			ushort LocalSeedVOffset = (count == 0) ? cks->first.first : cks->second.first;
			// argB (see KmerIntersectSR.hpp) == col == seqH
			ushort LocalSeedHOffset = (count == 0) ? cks->first.second : cks->second.second;

			// Get sequences
			std::string seqH;
			std::string seqV;

			seqan::assign(seqH, seqsh[i]);
			seqan::assign(seqV, seqsv[i]);

			uint lenH = seqH.length();
			uint lenV = seqV.length();

			// Get seed string
			std::string seedH = seqH.substr(LocalSeedHOffset, seed_length);
			std::string seedV = seqV.substr(LocalSeedVOffset, seed_length);

			std::string twinH = reversecomplement(seedH);

			if(twinH == seedV)
			{
				std::string twinseqH(seqH);

				std::reverse(std::begin(twinseqH), std::end(twinseqH));
				std::transform(std::begin(twinseqH), std::end(twinseqH), std::begin(twinseqH), complementbase);

				LocalSeedHOffset = twinseqH.length() - LocalSeedHOffset - seed_length;

				SeedInterface seed(LocalSeedHOffset, LocalSeedVOffset, seed_length); //LocalSeedHOffset + seed_length, LocalSeedVOffset + seed_length);

				// GGGG: here only accumulate stuff for the GPUs, don't perform alignment
				seeds.push_back(seed); // segfault origin might be around here 
				seqVs.push_back(seqV);
				seqHs.push_back(twinseqH);

				localRes.rc = true;
				xscores.push_back(localRes);
			}
			else
			{
				SeedInterface seed(LocalSeedHOffset, LocalSeedVOffset, seed_length); //LocalSeedHOffset + , LocalSeedVOffset + seed_length);

				// GGGG: here only accumulate stuff for the GPUs, don't perform alignment
				seeds.push_back(seed); // segfault origin might be around here (?)
				seqVs.push_back(seqV);
				seqHs.push_back(seqH);

				localRes.rc = false;
				xscores.push_back(localRes);
			}
		}

		auto end_time = std::chrono::system_clock::now();
    	add_time("XA:LoganPreprocess", (ms_t(end_time - start_time)).count());

		start_time = std::chrono::system_clock::now();

		// Call LOGAN only if noAlign is false
		if(!noAlign) 
		{ 
			// @GGGG-TODO: Check the parameter
			RunLoganAlign(seqHs, seqVs, seeds, xscores, xdrop, seed_length, myrank);
		}

		end_time = std::chrono::system_clock::now();
    	add_time("XA:LoganAlign", (ms_t(end_time - start_time)).count());

		start_time = std::chrono::system_clock::now();
		
		// @GGGG-TODO: this runs everything twice right now; pretty sure I can avoid this using CCS, I can add an option also run some benchmark
		// Compute stats
		if (count == 0)	// overwrite in the first seed
		{
			// @GGGG: keep the order for the post alignment evaluation (measure slowdown)
			// #pragma omp parallel for 
			for (uint64_t i = 0; i < npairs; ++i)
			{
				ai[i].xscore = xscores[i].score; 
				ai[i].rc     = xscores[i].rc;

                ai[i].begSeedH = xscores[i].begSeedH; 
                ai[i].endSeedH = xscores[i].endSeedH; 
                ai[i].begSeedV = xscores[i].begSeedV; 
                ai[i].endSeedV = xscores[i].endSeedV; 

				ai[i].seq_h_length = seqan::length(seqsh[i]);
				ai[i].seq_v_length = seqan::length(seqsv[i]);

				// this is a bit redundant since we can extract it from seed
				ai[i].seq_h_seed_length = ai[i].endSeedH - ai[i].begSeedH;
				ai[i].seq_v_seed_length = ai[i].endSeedV - ai[i].begSeedV;

				// GGGG: global idx over here to use in the FullDistVect for removing contained vertices/seqs
				ai[i].seq_h_g_idx = col_offset + std::get<1>(mattuples[lids[i]]);
    			ai[i].seq_v_g_idx = row_offset + std::get<0>(mattuples[lids[i]]);
			}
		}
		else
		{
			// @GGGG: keep the order for the post alignment evaluation (measure slowdown)
			// #pragma omp parallel for 
			for (uint64_t i = 0; i < npairs; ++i)
			{
				if (xscores[i].score > ai[i].xscore)
				{
					ai[i].xscore = xscores[i].score;
					ai[i].rc     = xscores[i].rc;

                    ai[i].begSeedH = xscores[i].begSeedH; 
                    ai[i].endSeedH = xscores[i].endSeedH; 
                    ai[i].begSeedV = xscores[i].begSeedV; 
                    ai[i].endSeedV = xscores[i].endSeedV; 

					// @GGGG: this is a bit redundant since we can extract it from seed
					ai[i].seq_h_seed_length = ai[i].endSeedH - ai[i].begSeedH;
					ai[i].seq_v_seed_length = ai[i].endSeedV - ai[i].begSeedV;
				}
			}
		}

		end_time = std::chrono::system_clock::now();
    	add_time("XA:ComputeStats", (ms_t(end_time - start_time)).count());
	}

	auto start_time = std::chrono::system_clock::now();
	std::vector<std::vector<int64_t>> ContainedSeqPerThread(numThreads);

	// Dump alignment info 
	// @GGGG: this should be fine parallel now (2/22 check here!)
	#pragma omp parallel
	{
	    #pragma omp for
		for (uint64_t i = 0; i < npairs; ++i)
		{
			// Only keep alignments that meet BELLA criteria
			bool passed = false;
			int tid = omp_get_thread_num();

			// GGGG: ai stores global idx to to store in ContainedSeqPerBatch
			// GGGG: in PostAlignDecision() we can mark as contained sequences as removable in ContainedSeqPerBatch and their local contained edges
			// GGGG: ContainedSeqPerBatch global indexes of contained sequences

			elba::CommonKmers *cks = std::get<2>(mattuples[lids[i]]);
			PostAlignDecision(ai[i], passed, ratioScoreOverlap, cks->overhang, cks->overhangT, cks->overlap, noAlign, ContainedSeqPerThread[tid]);

			if (passed)
			{
				// GGGG: store updated seed start/end position in the CommonKmers pairs (the semantics of these pairs change wrt the original semantics but that's okay)
				cks->first.first   = ai[i].begSeedV; 	// start on ver sequence
				cks->second.first  = ai[i].begSeedH;	    // start on hor sequence

				cks->first.second  = ai[i].endSeedV; 		// end on ver sequence
				cks->second.second = ai[i].endSeedH;		// end on hor sequence

				cks->lenv 	= ai[i].seq_v_length;
				cks->lenh 	= ai[i].seq_h_length;

				cks->score  = ai[i].xscore;
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
