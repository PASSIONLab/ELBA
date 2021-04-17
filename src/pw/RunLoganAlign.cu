/* Created by Giulia Guidi on 4/16/2021. */

#include "../../src/pw/GPULoganAligner.cpp"
#include "../../include/cuda/logan.cuh"

#define BATCH_SIZE 100000
#define MIN_OV_LEN 10000

void 
RunLoganAlign(vector<string>& seqHs, vector<string>& seqVs, vector<LSeed>& seeds, vector<loganResult>& xscores)
{
	ScoringSchemeL sscheme(1, -1, -1, -1);
	std::vector<ScoringSchemeL> scoring;
	scoring.push_back(sscheme);

	int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    omp_set_num_threads(deviceCount); // one OMP thread per GPU

	int AlignmentsToBePerformed = seeds.size();
	int numAlignmentsLocal = BATCH_SIZE * deviceCount; 

	//	Load balancer that divides the work in batches of 100K alignments
	for(int i = 0; i < AlignmentsToBePerformed; i += BATCH_SIZE * deviceCount)
	{
		if(AlignmentsToBePerformed < (i + BATCH_SIZE * deviceCount))
			numAlignmentsLocal = AlignmentsToBePerformed % (BATCH_SIZE * deviceCount);

		int* res = (int*)malloc(numAlignmentsLocal * sizeof(int));	

		std::vector<string>::const_iterator first_t = seqHs.begin() + i;
		std::vector<string>::const_iterator last_t  = seqHs.begin() + i + numAlignmentsLocal;
		std::vector<string> bseqHs(first_t, last_t);

		std::vector<string>::const_iterator first_q = seqVs.begin() + i;
		std::vector<string>::const_iterator last_q  = seqVs.begin() + i + numAlignmentsLocal;
		std::vector<string> bseqVs(first_q, last_q);

		std::vector<LSeed>::const_iterator first_s = seeds.begin() + i;
		std::vector<LSeed>::const_iterator last_s  = seeds.begin() + i + numAlignmentsLocal;
		std::vector<LSeed> bseeds(first_s, last_s);

		extendSeedL(bseeds, EXTEND_BOTHL, bseqHs, bseqVs, scoring, xdrop, seed_length, res, numAlignmentsLocal, deviceCount);

		for(int j = 0; j < numAlignmentsLocal; j++)
		{
			xscores[j+i].score = res[j];
			xscores[j+i].seed  = bseeds[j];
		}

		free(res);
	}
}