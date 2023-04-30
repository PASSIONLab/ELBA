#include "DistributedFastaData.hpp"
#include "Logger.hpp"

DistributedFastaData::DistributedFastaData(std::shared_ptr<FastaIndex> index) : index(index)
{
    Grid commgrid = index->getcommgrid();

    assert(commgrid->GetGridRows() == commgrid->GetGridCols());

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    int myrowid = commgrid->GetRankInProcCol();
    int mycolid = commgrid->GetRankInProcRow();
    int procdim = commgrid->GetGridRows();

    size_t numreads = index->gettotrecords();
    size_t readsperprocdim = numreads / procdim;

    isdiag = (myrowid == mycolid);

    rowstartid = myrowid * readsperprocdim;
    colstartid = mycolid * readsperprocdim;

    numrowreads = (myrowid == procdim-1)? (numreads - myrowid*readsperprocdim) : readsperprocdim;
    numcolreads = (mycolid == procdim-1)? (numreads - mycolid*readsperprocdim) : readsperprocdim;

    Logger logger(commgrid);
    logger() << "P(" << myrowid+1 << ", " << mycolid+1 << ") " << Logger::readrangestr(rowstartid, numrowreads) << "; " << Logger::readrangestr(colstartid, numcolreads);
    logger.Flush("DistributedFastaData::DistributedFastaData");

    auto requests = getremoterequests();
}

using FastaDataRequest = typename DistributedFastaData::FastaDataRequest;

void DistributedFastaData::getgridrequests(std::vector<FastaDataRequest>& myrequests, size_t globalstartid, size_t count, unsigned short rc) const
{
    Grid commgrid = index->getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    size_t totreads = index->gettotrecords();
    const auto& readdispls = index->getreaddispls();

    assert(count >= 1 && 0 <= globalstartid && globalstartid + count <= totreads);

    std::vector<int> nbrs;
    auto iditr = std::upper_bound(readdispls.cbegin(), readdispls.cend(), static_cast<MPI_Displ_type>(globalstartid));
    iditr--;

    assert(*iditr <= globalstartid);

    while (*iditr < globalstartid + count)
    {
        int sendrank = iditr - readdispls.cbegin();
        size_t reqstart = std::max(static_cast<size_t>(*iditr++), globalstartid);
        size_t reqend = std::min(static_cast<size_t>(*iditr), globalstartid + count);
        if (sendrank != myrank) myrequests.emplace_back(sendrank, myrank, reqstart, reqend-reqstart, rc);
    }
}

std::vector<FastaDataRequest> DistributedFastaData::getremoterequests() const
{
    std::vector<FastaDataRequest> allrequests, myrequests;

    getgridrequests(myrequests, rowstartid, numrowreads, 0);

    if (!isdiag) getgridrequests(myrequests, colstartid, numcolreads, 1);

    Logger logger(index->getcommgrid());

    logger() << "\n";
    for (auto itr = myrequests.begin(); itr != myrequests.end(); ++itr)
    {
        logger() << logger.rankstr(itr->requester) << " requests " << Logger::readrangestr(itr->offset, itr->count) << " from "
                 << logger.rankstr(itr->owner) << "; (rc=" << itr->rc << ")" << "\n";
    }
    logger.Flush("Remote requests:");
    return myrequests;
}
