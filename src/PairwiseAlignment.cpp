#include "PairwiseAlignment.hpp"
#include "DnaSeq.hpp"
#include "Logger.hpp"

CT<Overlap>::PSpParMat PairwiseAlignment(DistributedFastaData& dfd, CT<SharedSeeds>::PSpParMat& Bmat, int mat, int mis, int gap, int dropoff)
{
    auto index = dfd.getindex();
    auto commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();
    auto rowbuf = dfd.getrowbuf();
    auto colbuf = dfd.getcolbuf();

    size_t localnnzs = Bmat.seqptr()->getnnz();
    std::vector<CT<SharedSeeds>::ref_tuples> alignseeds;
    alignseeds.reserve(localnnzs);
    auto dcsc = Bmat.seqptr()->GetDCSC();

    uint64_t rowoffset = dfd.getrowstartid();
    uint64_t coloffset = dfd.getcolstartid();

    for (uint64_t i = 0; i < dcsc->nzc; ++i)
        for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
        {
            uint64_t localrow = dcsc->ir[j];
            uint64_t localcol = dcsc->jc[i];
            uint64_t globalrow = localrow + rowoffset;
            uint64_t globalcol = localcol + coloffset;

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

    std::vector<uint64_t> local_rowids, local_colids;
    std::vector<Overlap> overlaps;

    overlaps.reserve(nalignments);

    for (size_t i = 0; i < nalignments; ++i)
    {
        uint64_t localrow = std::get<0>(alignseeds[i]);
        uint64_t localcol = std::get<1>(alignseeds[i]);

        const DnaSeq& seqQ = (*rowbuf)[localrow];
        const DnaSeq& seqT = (*colbuf)[localcol];

        SeedPair len(static_cast<PosInRead>(seqQ.size()), static_cast<PosInRead>(seqT.size()));

        overlaps.emplace_back(len, std::get<2>(alignseeds[i])->getseeds()[0]);
        overlaps.back().extend_overlap(seqQ, seqT, mat, mis, gap, dropoff);

        local_rowids.push_back(localrow + rowoffset);
        local_colids.push_back(localcol + coloffset);
    }

    CT<uint64_t>::PDistVec drows(local_rowids, commgrid);
    CT<uint64_t>::PDistVec dcols(local_colids, commgrid);
    CT<Overlap>::PDistVec dvals(overlaps, commgrid);

    uint64_t numreads = index->gettotrecords();

    return CT<Overlap>::PSpParMat(numreads, numreads, drows, dcols, dvals, false);
}
