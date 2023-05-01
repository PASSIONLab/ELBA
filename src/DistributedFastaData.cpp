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

    size_t totreads = index->gettotrecords(); /* total reads in FASTA */
    const auto& readdispls = index->getreaddispls(); /* linear uniform read distribution displacements (across all processors) */

    assert(readdispls[nprocs] == totreads);

    /* assert that we are requesting more than zero reads, and that the range of reads exists within the FASTA */
    assert(count >= 1 && 0 <= globalstartid && globalstartid + count <= totreads);

    /*
     * We want to find the first processor rank that owns the read with id @globalstartid.
     * For each rank i in [0..nprocs-1], readdispls[i] is the number of reads owned across
     * processors [0..i), or the same: the global id of the first read owned by processor i.
     * readdispls[nprocs] is defined as the total number of reads in the FASTA.
     *
     * std:upper_bound(readdispls.cbegin(), readdispls.cend(), globalstartid) returns an
     * iterator pointing to the first element in the range [0..nprocs] such that
     * @globalstartid < element. This iterator is therefore pointing to the entry of
     * the first processor rank holding a read id greater than @globalstartid. It follows
     * that the processor rank just before it owns the read @globalstartid.
     */
    auto iditr = std::upper_bound(readdispls.cbegin(), readdispls.cend(), static_cast<MPI_Displ_type>(globalstartid));
    iditr--;

    /*
     * assert that the processor rank we think owns @globalstartid owns reads with
     * the same or lower id. If so, then we can be certain that we have found
     * the right starting rank because we already found from the std::upper_bound
     * computation that *(iditr+1) > @globalstartid.
     */
    assert(*iditr <= globalstartid);

    /*
     * @globalstartid + @count is the smallest index that we don't want for this
     * particular call to getgridrequests(). Therefore, we loop through the
     * processor ranks in the correct range and collect the information we
     * need from each one.
     */
    int requester = myrank;
    while (*iditr < globalstartid + count)
    {
        int owner = iditr - readdispls.cbegin();
        size_t reqstart = std::max(static_cast<size_t>(*iditr++), globalstartid);
        size_t reqend = std::min(static_cast<size_t>(*iditr), globalstartid + count);
        if (owner != requester) myrequests.emplace_back(owner, requester, reqstart, reqend-reqstart, rc);
        else myrequests.emplace_back(owner, requester, reqstart, reqend-reqstart, 3);
    }
}

/*
 * allrequests - allgathered requests
 * myrequests - requests originating from me
 */
void DistributedFastaData::getremoterequests(std::vector<FastaDataRequest>& allrequests, std::vector<FastaDataRequest>& myrequests) const
{
    Grid commgrid = index->getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();

    /*
     * Get row and column grid requests.
     */
    myrequests.resize(0);
    getgridrequests(myrequests, rowstartid, numrowreads, 0);
    if (!isdiag) getgridrequests(myrequests, colstartid, numcolreads, 1);

    Logger logger(commgrid);
    logger() << "\n";
    for (auto itr = myrequests.begin(); itr != myrequests.end(); ++itr)
        logger() << *itr << "\n";
    logger.Flush("myrequests");


    MPI_Count_type allrequestcount; /* total number of requests */
    MPI_Count_type myrequestcount = myrequests.size(); /* number of requests originating from me */

    std::vector<MPI_Count_type> requestcounts(nprocs); /* Allgatherv receive counts */
    std::vector<MPI_Displ_type> requestdispls(nprocs); /* Allgatherv receive displacements */

    /*
     * Globally collect the number of requests each processor wants to make.
     */
    requestcounts[myrank] = myrequestcount;
    MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_COUNT_TYPE, requestcounts.data(), 1, MPI_COUNT_TYPE, comm);

    /*
     * Compute allgatherv displacements.
     */
    std::exclusive_scan(requestcounts.begin(), requestcounts.end(), requestdispls.begin(), static_cast<MPI_Displ_type>(0));
    allrequestcount = requestdispls.back() + requestcounts.back();
    allrequests.resize(allrequestcount);

    /* Create MPI_Datatype for FastaDataRequest (needed for allgatherv) */
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

    /*
     * Globally collect all requests that are being made into @allrequests.
     */
    MPI_ALLGATHERV(myrequests.data(), myrequestcount, reqtype, allrequests.data(), requestcounts.data(), requestdispls.data(), reqtype, comm);
    MPI_Type_free(&reqtype);

    // std::sort(allrequests.begin(), allrequests.end(), [](const auto& a, const auto& b) { return a.owner < b.owner; });
}

void DistributedFastaData::blocking_read_exchange()
{
    std::vector<FastaDataRequest> allrequests, myrequests;
    getremoterequests(allrequests, myrequests);

    size_t totrowbases = 0;
    size_t totcolbases = 0;

    // rowbuf.reset(new DnaBuffer(totrowbases, getnumrowreads()));
    // colbuf.reset(new DnaBuffer(totcolbases, getnumcolreads()));
}
