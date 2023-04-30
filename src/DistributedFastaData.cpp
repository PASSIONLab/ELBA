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
        else myrequests.emplace_back(sendrank, myrank, reqstart, reqend-reqstart, 3);
    }
}

std::vector<FastaDataRequest> DistributedFastaData::getremoterequests() const
{
    Grid commgrid = index->getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();

    std::vector<FastaDataRequest> allrequests, myrequests;

    getgridrequests(myrequests, rowstartid, numrowreads, 0);
    if (!isdiag) getgridrequests(myrequests, colstartid, numcolreads, 1);

    MPI_Count_type allrequestcount, myrequestcount = myrequests.size();
    std::vector<MPI_Count_type> requestcounts(nprocs);
    std::vector<MPI_Displ_type> requestdispls(nprocs);

    requestcounts[myrank] = myrequestcount;
    MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_COUNT_TYPE, requestcounts.data(), 1, MPI_COUNT_TYPE, comm);

    std::exclusive_scan(requestcounts.begin(), requestcounts.end(), requestdispls.begin(), static_cast<MPI_Displ_type>(0));
    allrequestcount = requestdispls.back() + requestcounts.back();
    allrequests.resize(allrequestcount);

    MPI_Datatype reqtype;
    int blklens[5] = {1,1,1,1,1};
    MPI_Aint displs[5] = {offsetof(FastaDataRequest, owner),
                          offsetof(FastaDataRequest, requester),
                          offsetof(FastaDataRequest, offset),
                          offsetof(FastaDataRequest, count),
                          offsetof(FastaDataRequest, rc)};

    MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_SIZE_T, MPI_SIZE_T, MPI_UNSIGNED_SHORT};
    MPI_Type_create_struct(5, blklens, displs, types, &reqtype);
    MPI_Type_commit(&reqtype);
    MPI_ALLGATHERV(myrequests.data(), myrequestcount, reqtype, allrequests.data(), requestcounts.data(), requestdispls.data(), reqtype, comm);
    MPI_Type_free(&reqtype);

    std::sort(allrequests.begin(), allrequests.end(), [](const auto& a, const auto& b) { return a.owner < b.owner; });

    return allrequests;
}

void DistributedFastaData::blocking_read_exchange()
{

}
