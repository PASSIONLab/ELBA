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
    MPI_Type_commit(&reqtype);

    /*
     * Globally collect all requests that are being made into @allreqs;
     */
    MPI_ALLGATHERV(myreqs.data(), reqcounts[myrank], reqtype, allreqs.data(), reqcounts.data(), reqdispls.data(), reqtype, comm);
    MPI_Type_free(&reqtype);

    std::copy_if(allreqs.begin(), allreqs.end(), std::back_inserter(mysends), [&](const auto& req) { return req.owner == myrank; });
    mynumsends = mysends.size();

    std::vector<size_t> reqinfo(2*mynumreqs); /* even indices are number of reads, odd indices are buffer sizes */
    std::vector<MPI_Request> recvreqs(mynumreqs), sendreqs(mynumsends);

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        MPI_IRECV(reqinfo.data() + (2*i), 2, MPI_SIZE_T, myreqs[i].owner, 99, comm, recvreqs.data() + i);
    }

    logger() << "\n";
    for (size_t i = 0; i < mynumsends; ++i)
    {
        assert(mysends[i].owner == myrank);
        size_t offset = mysends[i].offset;
        size_t count = mysends[i].count;
        assert(index->getmyreaddispl() <= offset && offset + count <= index->getmyreaddispl() + index->getmyreadcount());
        size_t rangebufsize = mydna->getrangebufsize(offset - index->getmyreaddispl(), count);
        size_t sendbuf[2] = {count, rangebufsize};
        MPI_ISEND(sendbuf, 2, MPI_SIZE_T, mysends[i].requester, 99, comm, sendreqs.data() + i);
        logger() << "sent " << count << "/" << rangebufsize << " to " << logger.rankstr(mysends[i].requester) << " :: " << mysends[i] << "\n";
    }

    assert(mynumreqs <= std::numeric_limits<int>::max());
    assert(mynumsends <= std::numeric_limits<int>::max());

    MPI_Waitall(static_cast<int>(mynumsends), sendreqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(static_cast<int>(mynumreqs), recvreqs.data(), MPI_STATUSES_IGNORE);

    logger() << "\n";
    for (size_t i = 0; i < mynumreqs; ++i) logger() << "received " << reqinfo[2*i] << "/" << reqinfo[2*i + 1] << " from " << logger.rankstr(myreqs[i].owner) << " :: " << myreqs[i] << "\n";
    logger.Flush("Exchanges:");

    // std::vector<size_t> myreqcounts(mynumreqs), myreqcountsdispls(mynumreqs);
    // std::transform(myreqs.cbegin(), myreqs.cend(), myreqcounts.begin(), [](const auto& req) { return req.count; });
    // std::exclusive_scan(myreqcounts.cbegin(), myreqcounts.cend(), myreqcountsdispls.begin(), static_cast<size_t>(0));


    /* size_t totreqreadlens = std::accumulate(myreqs.begin(), myreqs.end(), static_cast<size_t>(0), [](size_t sum, const auto& req) { return sum + req.count; }); */



    return std::make_shared<DnaBuffer>(1);
}
