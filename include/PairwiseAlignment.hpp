#ifndef PAIRWISE_ALIGNMENT_H_
#define PAIRWISE_ALIGNMENT_H_

#include "common.h"
#include "DnaSeq.hpp"
#include "DistributedFastaData.hpp"
#include "SharedSeeds.hpp"
#include "Logger.hpp"

// void RunAlignments(DistributedFastaData& dfd, std::vector<CT<SharedSeeds>::ref_tuples>& alignin, std::vector<CT<Overlap>::ref_tuples>& alignout, int mat, int mis, int gap, int dropoff)
// {
    // auto rowbuf = dfd.getrowbuf();
    // auto colbuf = dfd.getcolbuf();

    // alignout.reserve(alignin.size());

    // auto commgrid = dfd.getindex()->getcommgrid();
    // Logger logger(commgrid);

    // for (auto itr = alignin.begin(); itr != alignin.end(); ++itr)
    // {
        // uint64_t localrow = std::get<0>(*itr);
        // uint64_t localcol = std::get<1>(*itr);
        // SharedSeeds *sseed = std::get<2>(*itr);

        // const DnaSeq& seqQ = (*rowbuf)[localrow];
        // const DnaSeq& seqT = (*colbuf)[localcol];

        // auto seedpairs = sseed->getseeds();
        // int numseeds = sseed->getnumstored();

        // for (int i = 0; i < numseeds; ++i)
        // {
            // alignout.emplace_back(seqQ, seqT, seedpairs[i]);
        // }
    // }

    // for (auto itr = alignout.begin(); itr != alignout.end(); ++itr)
    // {
        // itr->XDropSeedExtension(mat, mis, gap, dropoff);
    // }
// }

void PairwiseAlignment(DistributedFastaData& dfd, CT<SharedSeeds>::PSpParMat& Bmat, int mat, int mis, int gap, int dropoff)
{
    auto index = dfd.getindex();
    auto commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    size_t localnnzs = Bmat.seqptr()->getnnz();
    std::vector<CT<SharedSeeds>::ref_tuples> alignseeds;
    alignseeds.reserve(localnnzs);
    auto dcsc = Bmat.seqptr()->GetDCSC();

    for (uint64_t i = 0; i < dcsc->nzc; ++i)
        for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
        {
            uint64_t localrow = dcsc->ir[j];
            uint64_t localcol = dcsc->jc[i];
            uint64_t globalrow = localrow + dfd.getrowstartid();
            uint64_t globalcol = localcol + dfd.getcolstartid();

            if ((localrow < localcol) || (localrow <= localcol && globalrow < globalcol))
            {
                alignseeds.emplace_back(localrow, localcol, &(dcsc->numx[j]));
            }
        }

    size_t nalignments = alignseeds.size();
    size_t totalignments;

    MPI_ALLREDUCE(&nalignments, &totalignments, 1, MPI_SIZE_T, MPI_SUM, comm);

    Logger logger(commgrid);
    logger() << "performing " << nalignments << "/" << totalignments << " alignments";
    logger.Flush("Alignment Counts:");

    // std::vector<CT<Overlap>::ref_tuples> overlaps;

    // RunAlignments(dfd, alignments, overlaps, mat, mis, gap, dropoff);
}

#endif
