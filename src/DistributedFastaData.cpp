#include "DistributedFastaData.hpp"
#include "MPITypeHandler.hpp"
#include "Logger.hpp"
#include <limits>

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

    int requester = myrank; /* I'm the processor making a request */
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

    int owner = iditr - readdispls.cbegin();

    /*
     * assert that the processor rank we think owns @globalstartid owns reads with
     * the same or lower id. If so, then we can be certain that we have found
     * the right starting rank because we already found from the std::upper_bound
     * computation that *(iditr+1) > @globalstartid.
     */
    assert(readdispls[owner] <= globalstartid);

    /*
     * @globalstartid + @count is the smallest index that we don't want for this
     * particular call to getgridrequests(). Therefore, we loop through the
     * processor ranks in the correct range and collect the information we
     * need from each one.
     */
    while (readdispls[owner] < globalstartid + count)
    {
        assert(owner < nprocs);

         /* First read id stored on owning processor. */
        size_t ownerstartid = static_cast<size_t>(readdispls[owner]);

        /*
         * First readid stored on processor right after the owner,
         * or the total number of reads if the owner is the last
         * processor rank.
         */
        size_t nextstartid = static_cast<size_t>(readdispls[owner+1]);

         /* First readid we are requesting from owner. */
        size_t reqstart = std::max(ownerstartid, globalstartid);

         /* One past the last readid we are requesting from owner. */
        size_t reqend = std::min(nextstartid, globalstartid + count);

         /* If the owner and requester are the same, adjust the rcflag. */
        unsigned short rcflag = owner != requester? rc : rc + 2;
        myrequests.emplace_back(owner++, requester, reqstart, reqend - reqstart, rcflag);
    }
}

std::shared_ptr<DnaBuffer> DistributedFastaData::collect_row_sequences(std::shared_ptr<DnaBuffer> mydna)
{
    Grid commgrid = index->getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    Logger logger(commgrid);

    std::vector<FastaDataRequest> myreqs, mysends, allreqs;

    getgridrequests(myreqs, getrowstartid(), getnumrowreads(), 0);

    std::vector<MPI_Count_type> reqcounts(nprocs); /* Allgatherv receive counts */
    std::vector<MPI_Displ_type> reqdispls(nprocs); /* Allgatherv receive displacements */

    size_t totnumreqs, mynumsends;
    size_t mynumreqs = myreqs.size();

    /*
     * Globally collect the number of requests each processor wants to make.
     */
    reqcounts[myrank] = static_cast<MPI_Count_type>(mynumreqs);
    MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_COUNT_TYPE, reqcounts.data(), 1, MPI_COUNT_TYPE, comm);

    /*
     * Compute allgatherv displacements.
     */
    std::exclusive_scan(reqcounts.begin(), reqcounts.end(), reqdispls.begin(), static_cast<MPI_Displ_type>(0));
    totnumreqs = reqdispls.back() + reqcounts.back();
    allreqs.resize(totnumreqs);

    /*
     * Create MPI_Datatype for FastaDataRequest (needed for allgatherv).
     */
    MPI_Datatype reqtype;
    int blklens[5] = {1,1,1,1,1};
    MPI_Aint displs[5] = {offsetof(FastaDataRequest, owner), offsetof(FastaDataRequest, requester), offsetof(FastaDataRequest, offset), offsetof(FastaDataRequest, count), offsetof(FastaDataRequest, rc)};
    MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_SIZE_T, MPI_SIZE_T, MPI_UNSIGNED_SHORT};
    MPI_Type_create_struct(5, blklens, displs, types, &reqtype);
    MPITypeHandler rthandler(&reqtype);

    /*
     * Globally collect all requests that are being made into @allreqs;
     */
    MPI_ALLGATHERV(myreqs.data(), reqcounts[myrank], reqtype, allreqs.data(), reqcounts.data(), reqdispls.data(), reqtype, comm);

    std::copy_if(allreqs.begin(), allreqs.end(), std::back_inserter(mysends), [&](const auto& req) { return req.owner == myrank; });
    mynumsends = mysends.size();


    std::vector<size_t> reqbufsizes(mynumreqs);
    std::vector<MPI_Request> recvreqs(mynumreqs), sendreqs(mynumsends);

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        MPI_IRECV(reqbufsizes.data() + i, 1, MPI_SIZE_T, myreqs[i].owner, 99, comm, recvreqs.data() + i);
    }

    logger() << "\n";
    for (size_t i = 0; i < mynumsends; ++i)
    {
        assert(mysends[i].owner == myrank);
        size_t offset = mysends[i].offset;
        size_t count = mysends[i].count;
        assert(index->getmyreaddispl() <= offset && offset + count <= index->getmyreaddispl() + index->getmyreadcount());
        size_t rangebufsize = mydna->getrangebufsize(offset - index->getmyreaddispl(), count);
        MPI_ISEND(&rangebufsize, 1, MPI_SIZE_T, mysends[i].requester, 99, comm, sendreqs.data() + i);
        logger() << "sent " << rangebufsize << " to " << logger.rankstr(mysends[i].requester) << " :: " << mysends[i] << "\n";
    }

    assert(mynumreqs <= std::numeric_limits<int>::max());
    assert(mynumsends <= std::numeric_limits<int>::max());

    MPI_Waitall(static_cast<int>(mynumsends), sendreqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(static_cast<int>(mynumreqs), recvreqs.data(), MPI_STATUSES_IGNORE);

    logger() << "\n";
    for (size_t i = 0; i < mynumreqs; ++i) logger() << "received " << reqbufsizes[i] << " from " << logger.rankstr(myreqs[i].owner) << " :: " << myreqs[i] << "\n";
    logger.Flush("Exchanges:");

    return std::make_shared<DnaBuffer>(1);
}

/*
 * allrequests - allgathered requests
 * myrequests - requests originating from me
 */
//void DistributedFastaData::getremoterequests(std::vector<FastaDataRequest>& allrequests, std::vector<FastaDataRequest>& myrequests) const
//{
//    Grid commgrid = index->getcommgrid();
//    int nprocs = commgrid->GetSize();
//    int myrank = commgrid->GetRank();
//    MPI_Comm comm = commgrid->GetWorld();
//
//    /*
//     * Get row and column grid requests.
//     */
//    myrequests.resize(0);
//    getgridrequests(myrequests, rowstartid, numrowreads, 0);
//
//    MPI_Count_type allrequestcount; /* total number of requests */
//    MPI_Count_type myrequestcount = myrequests.size(); /* number of requests originating from me */
//
//    std::vector<MPI_Count_type> requestcounts(nprocs); /* Allgatherv receive counts */
//    std::vector<MPI_Displ_type> requestdispls(nprocs); /* Allgatherv receive displacements */
//
//    /*
//     * Globally collect the number of requests each processor wants to make.
//     */
//    requestcounts[myrank] = myrequestcount;
//    MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_COUNT_TYPE, requestcounts.data(), 1, MPI_COUNT_TYPE, comm);
//
//    /*
//     * Compute allgatherv displacements.
//     */
//    std::exclusive_scan(requestcounts.begin(), requestcounts.end(), requestdispls.begin(), static_cast<MPI_Displ_type>(0));
//    allrequestcount = requestdispls.back() + requestcounts.back();
//    allrequests.resize(allrequestcount);
//
//    /* Create MPI_Datatype for FastaDataRequest (needed for allgatherv) */
//    MPI_Datatype reqtype;
//    int blklens[5] = {1,1,1,1,1};
//    MPI_Aint displs[5] = {offsetof(FastaDataRequest, owner),
//                          offsetof(FastaDataRequest, requester),
//                          offsetof(FastaDataRequest, offset),
//                          offsetof(FastaDataRequest, count),
//                          offsetof(FastaDataRequest, rc)};
//
//    MPI_Datatype types[5] = {MPI_INT, MPI_INT, MPI_SIZE_T, MPI_SIZE_T, MPI_UNSIGNED_SHORT};
//    MPI_Type_create_struct(5, blklens, displs, types, &reqtype);
//    MPI_Type_commit(&reqtype);
//
//    /*
//     * Globally collect all requests that are being made into @allrequests.
//     */
//    MPI_ALLGATHERV(myrequests.data(), myrequestcount, reqtype, allrequests.data(), requestcounts.data(), requestdispls.data(), reqtype, comm);
//    MPI_Type_free(&reqtype);
//}


//void DistributedFastaData::blocking_read_exchange(std::shared_ptr<DnaBuffer> mydna)
//{
//    Grid commgrid = index->getcommgrid();
//    int nprocs = commgrid->GetSize();
//    int myrank = commgrid->GetRank();
//    MPI_Comm comm = commgrid->GetWorld();
//    Logger logger(commgrid);
//
//    std::vector<FastaDataRequest> allrequests, myrequests;
//    getremoterequests(allrequests, myrequests);
//
//    size_t numallrequests, nummyrequests, nummysends;
//    std::vector<size_t> recvbufsizes;
//    std::vector<MPI_Request> recvrequests, sendrequests;
//
//    nummyrequests = myrequests.size();
//    numallrequests = allrequests.size();
//    recvbufsizes.resize(nummyrequests);
//    recvrequests.resize(nummyrequests);
//
//    for (size_t i = 0; i < nummyrequests; ++i)
//    {
//        MPI_IRECV(recvbufsizes.data() + i, 1, MPI_SIZE_T, myrequests[i].owner, 99+myrequests[i].rc, comm, recvrequests.data() + i);
//    }
//
//    std::vector<std::tuple<size_t, int, unsigned short>> sendinfo;
//
//    for (size_t i = 0; i < numallrequests; ++i)
//    {
//        if (allrequests[i].owner == myrank)
//        {
//            size_t start = allrequests[i].offset - index->getmyreaddispl();
//            size_t sendbufsize = mydna->getrangebufsize(start, allrequests[i].count);
//            sendinfo.emplace_back(sendbufsize, allrequests[i].requester, allrequests[i].rc);
//        }
//    }
//
//    nummysends = sendinfo.size();
//    sendrequests.resize(nummysends);
//
//    for (size_t i = 0; i < nummysends; ++i)
//    {
//        size_t sendbufsize = std::get<0>(sendinfo[i]);
//        int dest = std::get<1>(sendinfo[i]);
//        unsigned short rc = std::get<2>(sendinfo[i]);
//        MPI_ISEND(&sendbufsize, 1, MPI_SIZE_T, dest, 99+rc, comm, sendrequests.data() + i);
//    }
//
//    assert(nummyrequests <= std::numeric_limits<int>::max());
//    assert(nummysends <= std::numeric_limits<int>::max());
//
//    MPI_Waitall(static_cast<int>(nummysends),    sendrequests.data(), MPI_STATUSES_IGNORE);
//    MPI_Waitall(static_cast<int>(nummyrequests), recvrequests.data(), MPI_STATUSES_IGNORE);
//
//    // rowbuf.reset(new DnaBuffer())
//}
