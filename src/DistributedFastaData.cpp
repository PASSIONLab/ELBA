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
    std::vector<MPI_Request> recvreqs(2*mynumreqs), sendreqs(mynumsends);

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        MPI_IRECV(reqinfo.data() + (2*i), 2, MPI_SIZE_T, myreqs[i].owner, 99, comm, recvreqs.data() + i);
    }

    std::vector<size_t> sendlens(mynumsends), sendbufsizes(mynumsends);
    std::vector<const uint8_t*> sendbufs(mynumsends);

    logger() << "\n";
    for (size_t i = 0; i < mynumsends; ++i)
    {
        assert(mysends[i].owner == myrank);
        size_t offset = mysends[i].offset;
        size_t localoffset = offset - index->getmyreaddispl();
        sendlens[i] = mysends[i].count;
        assert(index->getmyreaddispl() <= offset && offset + sendlens[i] <= index->getmyreaddispl() + index->getmyreadcount());
        sendbufsizes[i] = mydna->getrangebufsize(localoffset, sendlens[i]);
        sendbufs[i] = mydna->getbufoffset(localoffset);
        size_t sendbuf[2] = {sendlens[i], sendbufsizes[i]};
        MPI_ISEND(sendbuf, 2, MPI_SIZE_T, mysends[i].requester, 99, comm, sendreqs.data() + i);
        logger() << "sent " << sendlens[i] << "/" << sendbufsizes[i] << " to " << logger.rankstr(mysends[i].requester) << " :: " << mysends[i] << "\n";
    }

    assert(2*mynumreqs <= std::numeric_limits<int>::max());
    assert(2*mynumsends <= std::numeric_limits<int>::max());

    MPI_Waitall(static_cast<int>(mynumsends), sendreqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(static_cast<int>(mynumreqs), recvreqs.data(), MPI_STATUSES_IGNORE);

    logger() << "\n";
    for (size_t i = 0; i < mynumreqs; ++i) logger() << "received " << reqinfo[2*i] << "/" << reqinfo[2*i + 1] << " from " << logger.rankstr(myreqs[i].owner) << " :: " << myreqs[i] << "\n";
    logger.Flush("Exchanges:");

    std::vector<size_t> reqreadlendispls(mynumreqs+1);
    std::vector<size_t> reqbufdispls(mynumreqs+1);
    reqreadlendispls.front() = reqbufdispls.front() = 0;

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        reqreadlendispls[i+1] = reqreadlendispls[i] + reqinfo[2*i];
        reqbufdispls[i+1] = reqbufdispls[i] + reqinfo[2*i+1];
    }

    size_t reqnumreads = reqreadlendispls.back();
    size_t reqbufsize = reqbufdispls.back();
    size_t *reqreadlens = new size_t[reqnumreads];
    uint8_t *reqbuf = new uint8_t[reqbufsize];

    logger() << "\n";
    for (size_t i = 0; i < mynumreqs; ++i)
    {
        MPI_Count_type count = reqreadlendispls[i+1] - reqreadlendispls[i];
        MPI_Count_type bufsize = reqbufdispls[i+1] - reqbufdispls[i];

        logger() << "receiving " << count << " read lengths at displacement " << reqreadlendispls[i] << " from " << logger.rankstr(myreqs[i].owner) << "\n";
        logger() << "receiving " << bufsize << " buffer bytes at displacement " << reqbufdispls[i] << " from " << logger.rankstr(myreqs[i].owner) << "\n";
        MPI_IRECV(reqreadlens + reqreadlendispls[i], count, MPI_SIZE_T, myreqs[i].owner, 100, comm, recvreqs.data() + i);
        MPI_IRECV(reqbuf + reqbufdispls[i], bufsize, MPI_UINT8_T, myreqs[i].owner, 101, comm, recvreqs.data() + mynumreqs + i);
    }
    logger.Flush("posted received calls:");

    std::vector<size_t> myreadlens = index->getmyreadlens();

    logger() << "\n";
    for (size_t i = 0; i < mynumsends; ++i)
    {
        size_t localoffset = mysends[i].offset - index->getmyreaddispl();
        const size_t *sendreadlens = myreadlens.data() + localoffset;
        const uint8_t *sendbuf = mydna->getbufoffset(localoffset);
        logger() << "sending " << sendlens[i] << " read lengths to " << logger.rankstr(mysends[i].requester) << "\n";
        logger() << "sending " << sendbufsizes[i] << " buffer bytes to " << logger.rankstr(mysends[i].requester) << "\n";
        MPI_ISEND(sendreadlens, static_cast<MPI_Count_type>(sendlens[i]), MPI_SIZE_T, mysends[i].requester, 100, comm, sendreqs.data() + i);
        MPI_ISEND(sendbuf, static_cast<MPI_Count_type>(reqinfo[2*i+1]), MPI_UINT8_T, mysends[i].requester, 101, comm, sendreqs.data() + mynumsends + i);
    }
    logger.Flush("posted send calls:");

    std::vector<MPI_Status> sendstatuses(2*mynumsends), recvstatuses(2*mynumreqs);

    MPI_Waitall(static_cast<int>(2*mynumsends), sendreqs.data(), sendstatuses.data());

    logger() << "\n";
    for (size_t i = 0; i < 2*mynumsends; ++i)
    {
        logger() << "sendstatuses[" << i << "]" << sendstatuses[i].MPI_SOURCE << "/" << sendstatuses[i].MPI_TAG << "/" << sendstatuses[i].MPI_ERROR << "\n";
    }
    logger.Flush("Send Statuses:");

    MPI_Waitall(static_cast<int>(2*mynumreqs), recvreqs.data(), recvstatuses.data());
    logger() << "\n";
    for (size_t i = 0; i < 2*mynumreqs; ++i)
    {
        logger() << "recvstatuses[" << i << "]" << recvstatuses[i].MPI_SOURCE << "/" << recvstatuses[i].MPI_TAG << "/" << recvstatuses[i].MPI_ERROR << "\n";
    }
    logger.Flush("Receive Statuses:");

    return std::make_shared<DnaBuffer>(1);
}
