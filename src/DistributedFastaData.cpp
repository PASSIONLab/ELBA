#include "DistributedFastaData.hpp"
#include "Logger.hpp"
#include <limits>
#include <iomanip>

DistributedFastaData::DistributedFastaData(FastaIndex& index) : index(index)
{
    std::shared_ptr<CommGrid> commgrid = index.getcommgrid();

    assert(commgrid->GetGridRows() == commgrid->GetGridCols());

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    int myrowid = commgrid->GetRankInProcCol();
    int mycolid = commgrid->GetRankInProcRow();
    int procdim = commgrid->GetGridRows();

    size_t numreads = index.gettotrecords();
    size_t readsperprocdim = numreads / procdim;

    isdiag = (myrowid == mycolid);

    rowinfo.startid = myrowid * readsperprocdim;
    colinfo.startid = mycolid * readsperprocdim;

    rowinfo.numreads = (myrowid == procdim-1)? (numreads - myrowid*readsperprocdim) : readsperprocdim;
    colinfo.numreads = (mycolid == procdim-1)? (numreads - mycolid*readsperprocdim) : readsperprocdim;

    #if LOG_LEVEL >= 2
    Logger logger(commgrid);
    logger() << "P(" << myrowid+1 << ", " << mycolid+1 << ") " << Logger::readrangestr(rowinfo.startid, rowinfo.numreads) << "; " << Logger::readrangestr(colinfo.startid, colinfo.numreads);
    logger.Flush("DistributedFastaData::DistributedFastaData");
    #endif
}

using FastaDataRequest = typename DistributedFastaData::FastaDataRequest;

void DistributedFastaData::getgridrequests(std::vector<FastaDataRequest>& myrequests, size_t globalstartid, size_t count, unsigned short rc) const
{
    std::shared_ptr<CommGrid> commgrid = index.getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    Logger logger(commgrid);

    int requester = myrank; /* I'm the processor making a request */
    size_t totreads = index.gettotrecords(); /* total reads in FASTA */
    const auto& readdispls = index.getreaddispls(); /* linear uniform read distribution displacements (across all processors) */

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
        //unsigned short rcflag = owner != requester? rc : rc + 2;
        myrequests.emplace_back(owner++, requester, reqstart, reqend - reqstart, rc);
    }
}

void DistributedFastaData::collect_sequences(const DnaBuffer& mydna)
{
    collect_dim_sequences(mydna, rowinfo);
    collect_dim_sequences(mydna, colinfo);
}

void DistributedFastaData::collect_dim_sequences(const DnaBuffer& mydna, DimExchangeInfo& diminfo)
{
    assert(&diminfo == &rowinfo || &diminfo == &colinfo);

    unsigned short rc = &diminfo == &rowinfo? 0 : 1;
    std::shared_ptr<CommGrid> commgrid = index.getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    Logger logger(commgrid);

    std::vector<FastaDataRequest> myreqs, mysends, allreqs;

    getgridrequests(myreqs, diminfo.startid, diminfo.numreads, rc);

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
    diminfo.recvreqs.resize(2*mynumreqs);

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        MPI_IRECV(reqinfo.data() + (2*i), 2, MPI_SIZE_T, myreqs[i].owner, 100+myreqs[i].rc, comm, diminfo.recvreqs.data() + i);
    }

    std::vector<size_t> sendlens(mynumsends), sendbufsizes(mynumsends);
    std::vector<const uint8_t*> sendbufs(mynumsends);

    for (size_t i = 0; i < mynumsends; ++i)
    {
        assert(mysends[i].owner == myrank);
        size_t offset = mysends[i].offset;
        size_t localoffset = offset - index.getmyreaddispl();
        sendlens[i] = mysends[i].count;
        assert(index.getmyreaddispl() <= offset && offset + sendlens[i] <= index.getmyreaddispl() + index.getmyreadcount());
        sendbufsizes[i] = mydna.getrangebufsize(localoffset, sendlens[i]);
        sendbufs[i] = mydna.getbufoffset(localoffset);
        size_t sendbuf[2] = {sendlens[i], sendbufsizes[i]};
        MPI_SEND(sendbuf, 2, MPI_SIZE_T, mysends[i].requester, 100+mysends[i].rc, comm);
    }

    assert(2*mynumreqs <= std::numeric_limits<int>::max());
    assert(2*mynumsends <= std::numeric_limits<int>::max());

    /*
     * Note that we don't use 2*mynum{sends,reqs} below for the count.
     */
    MPI_Waitall(static_cast<int>(mynumreqs), diminfo.recvreqs.data(), MPI_STATUSES_IGNORE);

    std::vector<size_t> reqreadlendispls(mynumreqs+1);
    std::vector<size_t> reqbufdispls(mynumreqs+1);
    reqreadlendispls.front() = reqbufdispls.front() = 0;

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        reqreadlendispls[i+1] = reqreadlendispls[i] + reqinfo[2*i];
        reqbufdispls[i+1] = reqbufdispls[i] + reqinfo[2*i+1];
    }

    diminfo.reqnumreads = reqreadlendispls.back();
    diminfo.reqbufsize = reqbufdispls.back();
    diminfo.reqreadlens.reset(new size_t[diminfo.reqnumreads]);
    diminfo.reqbuf.reset(new uint8_t[diminfo.reqbufsize]);

    for (size_t i = 0; i < mynumreqs; ++i)
    {
        MPI_Count_type count = reqreadlendispls[i+1] - reqreadlendispls[i];
        MPI_Count_type bufsize = reqbufdispls[i+1] - reqbufdispls[i];
        MPI_IRECV(diminfo.reqreadlens.get() + reqreadlendispls[i], count, MPI_SIZE_T, myreqs[i].owner, 200+myreqs[i].rc, comm, diminfo.recvreqs.data() + i);
        MPI_IRECV(diminfo.reqbuf.get() + reqbufdispls[i], bufsize, MPI_UINT8_T, myreqs[i].owner, 300+myreqs[i].rc, comm, diminfo.recvreqs.data() + mynumreqs + i);
    }

    std::vector<size_t> myreadlens = index.getmyreadlens();

    for (size_t i = 0; i < mynumsends; ++i)
    {
        size_t localoffset = mysends[i].offset - index.getmyreaddispl();
        const size_t *sendreadlens = myreadlens.data() + localoffset;
        const uint8_t *sendbuf = mydna.getbufoffset(localoffset);
        MPI_SEND(sendreadlens, static_cast<MPI_Count_type>(sendlens[i]), MPI_SIZE_T, mysends[i].requester, 200+mysends[i].rc, comm);
        MPI_SEND(sendbuf, static_cast<MPI_Count_type>(sendbufsizes[i]), MPI_UINT8_T, mysends[i].requester, 300+mysends[i].rc, comm);
    }
}

void DistributedFastaData::wait()
{
    MPI_Waitall(static_cast<int>(rowinfo.recvreqs.size()), rowinfo.recvreqs.data(), MPI_STATUSES_IGNORE);
    MPI_Waitall(static_cast<int>(colinfo.recvreqs.size()), colinfo.recvreqs.data(), MPI_STATUSES_IGNORE);

    rowbuf.reset(new DnaBuffer(rowinfo.reqbufsize, rowinfo.reqnumreads, rowinfo.reqbuf.release(), rowinfo.reqreadlens.get()));
    colbuf.reset(new DnaBuffer(colinfo.reqbufsize, colinfo.reqnumreads, colinfo.reqbuf.release(), colinfo.reqreadlens.get()));
}

std::string getgridfname(char const *fname_prefix, int rank, bool rc)
{
    std::ostringstream ss;
    ss << fname_prefix << (rc? "col" : "row") << rank+1 << ".txt";
    return ss.str();
}

void DistributedFastaData::write_grid_sequences(char const *fname_prefix) const
{
    auto commgrid = index.getcommgrid();
    MPI_Comm comm = commgrid->GetWorld();
    MPI_Comm rowcomm = commgrid->GetRowWorld();
    MPI_Comm colcomm = commgrid->GetColWorld();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    int myrowid = commgrid->GetRankInProcCol(); /* same as MPI_Comm_rank(colcomm, &myrowid) */
    int mycolid = commgrid->GetRankInProcRow(); /* same as MPI_Comm_rank(rowcomm, &mycolid) */
    int procdim = commgrid->GetGridRows();

    std::string myrowfname = getgridfname(fname_prefix, myrowid, 0);
    std::string mycolcontents = colbuf->getasciifilecontents();
    MPI_Count mycolcount = mycolcontents.size();
    MPI_File fh_rows;
    MPI_File_open(rowcomm, myrowfname.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_rows);
    MPI_File_write_ordered(fh_rows, mycolcontents.c_str(), mycolcount, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh_rows);

    std::string mycolfname = getgridfname(fname_prefix, mycolid, 1);
    std::string myrowcontents = rowbuf->getasciifilecontents();
    MPI_Count myrowcount = myrowcontents.size();
    MPI_File fh_cols;
    MPI_File_open(colcomm, mycolfname.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh_cols);
    MPI_File_write_ordered(fh_cols, myrowcontents.c_str(), myrowcount, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh_cols);

    MPI_Barrier(comm);
}

