#ifndef PAIRWISE_ALIGNMENT_H_
#define PAIRWISE_ALIGNMENT_H_

#include "common.h"
#include "DnaSeq.hpp"
#include "DistributedFastaData.hpp"
#include "SharedSeeds.hpp"
#include "Logger.hpp"

void RunAlignments(DistributedFastaData& dfd, std::vector<CT<SharedSeeds>::ref_tuples>& localtuples, int mat, int mis, int gap, int dropoff)
{
    auto rowbuf = dfd.getrowbuf();
    auto colbuf = dfd.getcolbuf();

    auto commgrid = dfd.getindex()->getcommgrid();
    Logger logger(commgrid);

    logger() << "\n";
    for (auto itr = localtuples.begin(); itr != localtuples.end(); ++itr)
    {
        uint64_t localrow = std::get<0>(*itr);
        uint64_t localcol = std::get<1>(*itr);
        SharedSeeds *sseed = std::get<2>(*itr);

        int rowseqlen = (*rowbuf)[localrow].size();
        int colseqlen = (*colbuf)[localcol].size();

        auto seedpairs = sseed->getseeds();
        int numseeds = sseed->getnumstored();

        for (int i = 0; i < numseeds; ++i)
        {
            const SeedPair& p = seedpairs[i];
            logger() << localrow << "\t" << localcol << "\t" << std::get<0>(p) << "\t" << std::get<1>(p) << "\t" << rowseqlen << "\t" << colseqlen << "\n";
        }
    }
    logger.Flush("Alignments to perform:");
}

void PairwiseAlignment(DistributedFastaData& dfd, CT<SharedSeeds>::PSpParMat& Bmat, int mat, int mis, int gap, int dropoff)
{
    auto index = dfd.getindex();
    auto commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    size_t localnnzs = Bmat.seqptr()->getnnz();
    std::vector<CT<SharedSeeds>::ref_tuples> localtuples, alignments;
    localtuples.reserve(localnnzs);
    auto dcsc = Bmat.seqptr()->GetDCSC();

    for (uint64_t i = 0; i < dcsc->nzc; ++i)
        for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
            localtuples.emplace_back(dcsc->ir[j], dcsc->jc[i], &dcsc->numx[j]);

    size_t nalignments = 0;

    for (size_t i = 0; i < localnnzs; ++i)
    {
        uint64_t localrow = std::get<0>(localtuples[i]);
        uint64_t localcol = std::get<1>(localtuples[i]);
        SharedSeeds *sseed = std::get<2>(localtuples[i]);

        uint64_t globalrow = localrow + dfd.getrowstartid();
        uint64_t globalcol = localcol + dfd.getcolstartid();

        if ((localrow < localcol) || (localrow <= localcol && globalrow < globalcol))
        {
            nalignments++;
            alignments.emplace_back(localrow, localcol, sseed);
        }
    }

    size_t totalignments;
    MPI_ALLREDUCE(&nalignments, &totalignments, 1, MPI_SIZE_T, MPI_SUM, comm);

    Logger logger(commgrid);
    logger() << "performing " << nalignments << "/" << totalignments << " alignments";
    logger.Flush("Alignment Counts:");

    RunAlignments(dfd, alignments, mat, mis, gap, dropoff);
}

#endif
