/* Created by Giulia Guidi on 4/20/2021. */

#include "RunLoganAligner.hpp"
#include "logan.hpp"

using namespace std;

void 
RunLoganAlign(vector<string>& seqHs, vector<string>& seqVs, 
	vector<SeedInterface>& SeedInterfaceSet, vector<LoganResult>& xscores, int& xdrop, ushort& seed_length)
{
	ScoringSchemeL sscheme(1, -1, -1, -1);
	std::vector<ScoringSchemeL> scoring;
	scoring.push_back(sscheme);

	int deviceCount;
	cudaGetDeviceCount(&deviceCount); // 1 MPI process to many GPUs 

	std::cout << deviceCount << " GPUs" << std::endl;

	int AlignmentsToBePerformed = SeedInterfaceSet.size();
	int numAlignmentsLocal = BATCH_SIZE * deviceCount; 

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

		std::vector<SeedInterface>::const_iterator first_s = SeedInterfaceSet.begin() + i;
		std::vector<SeedInterface>::const_iterator last_s  = SeedInterfaceSet.begin() + i + numAlignmentsLocal;
		
		std::vector<SeedInterface> bSeedInterface(first_s, last_s);

		std::vector<LSeed> bLSeedSet;
		for(int k = 0; k < (int)bSeedInterface.size(); k++)
		{
			LSeed lseed(bSeedInterface[k]);
			bLSeedSet.push_back(lseed); // segfault origin might be around here 
		}

		extendSeedL(bLSeedSet, EXTEND_BOTHL, bseqHs, bseqVs, scoring, xdrop, seed_length, res, numAlignmentsLocal, deviceCount, GPU_THREADS);

		for(int j = 0; j < numAlignmentsLocal; j++)
		{
			xscores[j+i].score     = res[j];

			xscores[j+i].begSeedH  = getBeginPositionH(bLSeedSet[j]);
			xscores[j+i].begSeedV  = getBeginPositionV(bLSeedSet[j]);

			xscores[j+i].endSeedH  = getEndPositionH(bLSeedSet[j]);
			xscores[j+i].endSeedV  = getEndPositionV(bLSeedSet[j]);
        }

		free(res);
	}
}



