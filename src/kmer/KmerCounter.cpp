//
// Originally authored and developed for short read analysis by: Aydin Buluc, Evangelos Georganas, and Rob Egan.
// Extended for long read analysis by Marquita Ellis and Giulia Guidi.
//

#include <numeric>
#include "KmerOps.cpp"
#include "../../include/NearestKmers2.hpp"

using namespace std;
using namespace dibella;

/////////////////////////////////////////////
// StoreReadName                           //
///////////////////////////////////////////// 

inline void StoreReadName(const ReadId& readIndex, std::string name, std::unordered_map<ReadId, std::string>& readNameMap)
{
  ASSERT(readNameMap.count(readIndex) == 0, "Rank "+ to_string(MYTHREAD) + ": collision in readNameMap on key = " + to_string(readIndex) + ", count is " + to_string(readNameMap.count(readIndex)));
  readNameMap[readIndex] = name;
}

/////////////////////////////////////////////
// countTotalKmersAndCleanHash             //
///////////////////////////////////////////// 

void countTotalKmersAndCleanHash()
{
    int64_t hashsize = 0;
    int64_t maxcount = 0;
    int64_t globalmaxcount = 0;

    /*! GGGG: where is kmercounts being filled? */
    for(auto itr = kmercounts->begin(); itr != kmercounts->end(); ++itr)
    {

      int allcount = get<2>(itr->second);
      if(allcount > maxcount)
        maxcount = allcount;

      nonerrorkmers += allcount;
      ++hashsize;
    }

    LOGF("my hashsize = %lld, nonerrorkmers = %lld, maxcount = %lld\n", (lld) hashsize, (lld) nonerrorkmers, (lld) maxcount);

    int64_t totalnonerror;
    int64_t distinctnonerror;

    CHECK_MPI(MPI_Reduce(&nonerrorkmers, &totalnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    CHECK_MPI(MPI_Reduce(&hashsize,   &distinctnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));
    CHECK_MPI(MPI_Allreduce(&maxcount,&globalmaxcount,   1, MPI_LONG_LONG, MPI_MAX,    MPI_COMM_WORLD));

    if(myrank == 0)
    {
        cout << "Counting finished " << endl;
        cout << __FUNCTION__ << ": Kmerscount hash includes " << distinctnonerror << " distinct elements" << endl;
        cout << __FUNCTION__ << ": Kmerscount non error kmers count is " << totalnonerror << endl;
        cout << __FUNCTION__ << ": Global max count is " << globalmaxcount << endl;
        cout << __FUNCTION__ << ": Large count histogram is of size " << HIGH_NUM_BINS << endl;

        /*! GGGG: define DIAG when code is clean and working */
        // ADD_DIAG("%lld", "distinct_non_error_kmers", (lld) distinctnonerror);
        // ADD_DIAG("%lld", "total_non_error_kmers", (lld) totalnonerror);
        // ADD_DIAG("%lld", "global_max_count", (lld) globalmaxcount);
    }

    /*! GGGG: heavy hitters part removed for now */  
    if (globalmaxcount == 0)
    {
        /*! GGGG: error message and terminate when code is clean and working */
        // SDIE("There were no kmers found, perhaps your KLEN (%d) is longer than your reads?", KLEN);
    }

    /* Reset */
    nonerrorkmers = 0;
    distinctnonerror = 0;

    int64_t overonecount = 0;

    auto itr = kmercounts->begin();
    while(itr != kmercounts->end())
    {
        int allcount =  get<2>(itr->second);
        if(allcount < ERR_THRESHOLD || (maxKmerFreq > 0 && allcount > maxKmerFreq))
        {
            --hashsize;
            itr = kmercounts->erase(itr);
        }
        else
        {
            nonerrorkmers += allcount;
            distinctnonerror++;
            ++itr;
        }

        if(allcount > 1)
        {
            overonecount += allcount;
        }
    }

    CHECK_MPI(MPI_Reduce(&nonerrorkmers, &totalnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD)); 
    CHECK_MPI(MPI_Reduce(&hashsize, &distinctnonerror, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));

    if(myrank == 0)
    {
        cout << __FUNCTION__ << ": Erroneous count < " << ERR_THRESHOLD  << " and high frequency > "<< maxKmerFreq <<" cases removed " << endl;
        cout << __FUNCTION__ << ": Kmerscount hash includes " << distinctnonerror << " distinct elements" << endl;
        cout << __FUNCTION__ << ": Kmerscount non error kmers count is " << totalnonerror << endl;

        // ADD_DIAG("%lld", "distinct_non_error_kmers", (lld) distinctnonerror);
        // ADD_DIAG("%lld", "total_non_error_kmers", (lld) totalnonerror);
        // ADD_DIAG("%lld", "global_max_count", (lld) globalmaxcount);
    }
}

/////////////////////////////////////////////
// DealWithInMemoryData                    //
/////////////////////////////////////////////

/*! GGGG: kmercounts (global variable) filled in DealWithInMemoryData */
/*  At this point, no kmers include anything other than uppercase 'A/C/G/T' */
void DealWithInMemoryData(VectorKmer& mykmers, int pass, struct bloom* bm, VectorReadId myreadids, VectorPos mypositions)
{
	LOGF("Dealing with memory pass %d: mykmers=%lld\n", pass, (lld) mykmers.size());

	if (pass == 2)
    {
		ASSERT(myreadids.size() == mykmers.size(), "");
		ASSERT(mypositions.size() == mykmers.size(), "");
    }
    else
    {
        ASSERT(pass == 1, "");
        ASSERT(myreadids.size() == 0, "");
        ASSERT(mypositions.size() == 0, "");
    }

    if (pass == 1)
    {
        assert(bm);
        size_t count = mykmers.size();
        for (size_t i = 0; i < count; ++i)
        {
            /* There will be a second pass, just insert into bloom, and map, but do not count */
            // KmerInfo ki(mykmers[i]);
            // ki.checkBloomAndInsert(bm, false);
        }
    }
    else
    {
        ASSERT(pass == 2, "");
        size_t count = mykmers.size();
        for(size_t i=0; i < count; ++i)
        {
            // KmerInfo ki(mykmers[i], myreadids[i], mypositions[i]);
			// ASSERT(!bm, "");
			// ki.includeCount(true);
        }
    }
}

/////////////////////////////////////////////
// ExchangePass                            //
/////////////////////////////////////////////

double Exchange(VectorVectorKmer& outgoing, VectorVectorReadId& readids, VectorVectorPos& positions, VectorVectorChar& extreads,
              VectorKmer& mykmers, VectorReadId& myreadids, VectorPos& mypositions, int pass, Buffer scratch1, Buffer scratch2)
{
    double totexch = MPI_Wtime();
    double perftime = 0.0;

    /*
     * Count and exchange number of bytes being sent
     * First pass: just k-mer (instances)
     * Second pass: each k-mer (instance) with its source read (ID) and position
     */

    size_t bytesperkmer  = Kmer::numBytes();
    size_t bytesperentry = bytesperkmer + (pass == 2 ? sizeof(ReadId) + sizeof(PosInRead) : 0);

    int* sendcnt = new int[nprocs];

    for(int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = (int) outgoing[i].size() * bytesperentry;

        // if (pass == 2)
        // {
        //     ASSERT( outgoing[i].size() == readids[i].size(), "" );
        //     ASSERT( outgoing[i].size() == positions[i].size(), "" );
        // }
        // else
        // {
        //     ASSERT (readids[i].size() == 0, "");
        //     ASSERT (positions[i].size() == 0, "");
        // }
    }

    int* sdispls = new int[nprocs];
    int* rdispls = new int[nprocs];
    int* recvcnt = new int[nprocs];

    /* Share the request counts */
    CHECK_MPI(MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD));  

    sdispls[0] = 0;
    rdispls[0] = 0;

    for(int i=0; i<nprocs-1; ++i)
    {
        if (sendcnt[i] < 0 || recvcnt[i] < 0)
        {
            cerr << myrank << " detected overflow in Alltoall" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        sdispls[i+1] = sdispls[i] + sendcnt[i];
        rdispls[i+1] = rdispls[i] + recvcnt[i];

        if (sdispls[i + 1] < 0 || rdispls[i + 1] < 0)
        {
            cerr << myrank << " detected overflow in Alltoall" << endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    int64_t totsend = accumulate(sendcnt, sendcnt+nprocs, static_cast<int64_t>(0));
    if (totsend < 0)
    {
        cerr << myrank << " detected overflow in totsend calculation, line" << __LINE__ << endl;
    }

    int64_t totrecv = accumulate(recvcnt, recvcnt+nprocs, static_cast<int64_t>(0));
    if (totrecv < 0)
    {
        cerr << myrank << " detected overflow in totrecv calculation, line" << __LINE__ << endl;
    }

    DBG("totsend = %lld totrecv = %lld\n", (lld) totsend, (lld) totrecv);

    /* It's gonna exit if totsend is negative */
    growBuffer(scratch1, sizeof(uint8_t) * totsend); 
    uint8_t * sendbuf = (uint8_t*) getStartBuffer(scratch1);

    for(int i = 0; i < nprocs; ++i)
    {
        size_t nkmers2send   = outgoing[i].size();
        uint8_t * addrs2fill = sendbuf+sdispls[i];

        for(size_t j = 0; j < nkmers2send; ++j)
        {
            ASSERT(addrs2fill == sendbuf+sdispls[i] + j*bytesperentry,"");
            (outgoing[i][j]).copyDataInto(addrs2fill);

            if (pass == 2)
            {
                ReadId* ptrRead = (ReadId*) (addrs2fill + bytesperkmer);
                ptrRead[0] = readids[i][j];
                PosInRead* ptrPos = (PosInRead*) (addrs2fill + bytesperkmer + sizeof(ReadId));
                ptrPos[0] = positions[i][j];
            }
            addrs2fill += bytesperentry;
        }

        outgoing[i].clear();
        readids[i].clear();
        positions[i].clear();
        extquals[i].clear();
        extreads[i].clear();
    }

    growBuffer(scratch2, sizeof(uint8_t) * totrecv);
    uint8_t * recvbuf = (uint8_t*) getStartBuffer(scratch2);

    double texch = 0.0 - MPI_Wtime();
    CHECK_MPI(MPI_Alltoallv(sendbuf, sendcnt, sdispls, MPI_BYTE, recvbuf, recvcnt, rdispls, MPI_BYTE, MPI_COMM_WORLD));
    texch += MPI_Wtime();

    /******* Performance report *******/
    perftime = MPI_Wtime();

    const int SND = 0;
    const int RCV = 1;

    int64_t local_counts[2];

    local_counts[SND] = totsend;
    local_counts[RCV] = totrecv;

    int64_t global_mins[2] = {0,0};
    CHECK_MPI(MPI_Reduce(&local_counts, &global_mins, 2, MPI_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD));

    int64_t global_maxs[2] = {0,0};
    CHECK_MPI(MPI_Reduce(&local_counts, &global_maxs, 2, MPI_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD));

    double global_min_time = 0.0;
    CHECK_MPI(MPI_Reduce(&texch, &global_min_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD));

    double global_max_time = 0.0;
    CHECK_MPI(MPI_Reduce(&texch, &global_max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));

    serial_printf("KmerMatch:%s exchange iteration %d pass %d: sent min %lld bytes, sent max %lld bytes, recv min %lld bytes, recv max %lld bytes, in min %.3f s, max %.3f s\n",
        __FUNCTION__, iter, pass, global_mins[SND], global_maxs[SND], global_mins[RCV], global_maxs[RCV], global_min_time, global_max_time);

    perftime = MPI_Wtime()-perftime;
    /*************************************/

    uint64_t nkmersrecvd = totrecv / bytesperentry;

    for(uint64_t i = 0; i < nkmersrecvd; ++i) 
    {
        Kmer kk;
        kk.copyDataFrom(recvbuf + (i * bytesperentry)); 
        mykmers.push_back(kk);

        if (pass == 2)
        {
            ReadId *ptr = (ReadId*) (recvbuf + (i * bytesperentry) + bytesperkmer);
            ASSERT(ptr[0] > 0, "");
            myreadids.push_back(ptr[0]);
            PosInRead *posPtr = (PosInRead*) (recvbuf + (i * bytesperentry) + bytesperkmer + sizeof(ReadId));
            mypositions.push_back(posPtr[0]);
        }
    }

    DBG("DeleteAll: recvcount = %lld, sendct = %lld\n", (lld) recvcnt, (lld) sendcnt);
    DeleteAll(rdispls, sdispls, recvcnt, sendcnt);

    iter++;
    totexch = MPI_Wtime() - totexch - perftime;

    // MPI_Pcontrol(-1,"Exchange");
    return totexch;
}

/////////////////////////////////////////////
// FinishPackPass1                         //
/////////////////////////////////////////////

/* The bloom filter pass; extensions are ignored */
inline size_t FinishPackPass1(VectorVectorKmer& outgoing, Kmer& kmerreal)
{
    /* whichever one is the representative */
    uint64_t myhash = kmerreal.hash(); 

    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());
    
    outgoing[owner].push_back(kmerreal);

    return outgoing[owner].size();
}

/////////////////////////////////////////////
// FinishPackPass2                         //
/////////////////////////////////////////////

/* The hash table pass; extensions are important */
inline size_t FinishPackPass2(VectorVectorKmer& outgoing, VectorVectorReadId& readids, VectorVectorPos& positions, 
    VectorVectorChar& extreads, Kmer& kmerreal, ReadId readId, PosInRead pos)
{
    assert(kmerreal == kmerreal.rep());

    /* whichever one is the representative */
    uint64_t myhash = kmerreal.hash();

    double range = static_cast<double>(myhash) * static_cast<double>(nprocs);
    size_t owner = range / static_cast<double>(numeric_limits<uint64_t>::max());

    size_t location = 0, maxsize = 0;

    /*! GGGG: find this information */
    /* Count here */
    outgoing[owner].push_back(kmerreal);
    readids[owner].push_back(readId);
    positions[owner].push_back(pos);

    return outgoing[owner].size();
}

/////////////////////////////////////////////
// PackEndsKmer function                   //
/////////////////////////////////////////////

size_t PackEndsKmer(string& seq, int j, Kmer& kmerreal, ReadId readid, PosInRead pos, VectorVectorKmer& outgoing,
        VectorVectorReadId& readids, VectorVectorPos& positions, VectorVectorChar& extreads, 
        int pass, int lastCountedBase, int KLEN)
{
    bool isCounted = lastCountedBase >= j + KLEN;
    size_t procSendCount;

    assert(seq.size() >= j + KLEN);
    assert(seq.substr(j, KLEN).find('N') == std::string::npos);
    assert(kmerreal == Kmer(seq.c_str() + j));

    if (pass == 1)
    {
        if (!isCounted) return 0;
        kmerreal = kmerreal.rep();
        procSendCount = FinishPackPass1(outgoing, kmerreal);
    }
    /* Otherwise we don't care about the extensions */
    else if (pass == 2)   
    {
        Kmer kmertwin = kmerreal.twin();
        /* The real k-mer is not lexicographically smaller */
        if(kmertwin < kmerreal) 
        {
            kmerreal = kmertwin;
        }
        procSendCount = FinishPackPass2(outgoing, readids, positions, extreads, kmerreal, readid, pos);
    }
    return procSendCount;
}

/////////////////////////////////////////////
// ParseNPack                              //
/////////////////////////////////////////////

/* Kmer is of length k */
/* HyperLogLog counting, bloom filtering, and std::maps use Kmer as their key */
size_t ParseNPack(FastaData& lfd, VectorVectorKmer& outgoing, VectorVectorReadId& readids, VectorVectorPos& positions, 
    ReadId& startReadIndex, VectorVectorChar& extreads, std::unordered_map<ReadId, std::string>& readNameMap, int pass, size_t offset)
{
    int len, start_offset, end_offset_inclusive;

    /* offset left over from orginal piece of codes */
    size_t nreads = lfd->local_count();
    size_t maxsending = 0, kmersthisbatch = 0, nskipped = 0;
    size_t bytesperkmer  = Kmer::numBytes();
    size_t bytesperentry = bytesperkmer + 4;
    size_t memthreshold  = (MAX_ALLTOALL_MEM/nprocs) * 2;

    /*! GGGG: where this startReadIndex come from? */
    ReadId readIndex = startReadIndex;

    /*! GGGG: make sure name are consistent with ids */
    for (uint64_t lreadidx = offset; lreadidx < nreads; ++lreadidx)
    {
        buff = lfd->get_sequence_id(lreadidx, len, start_offset, end_offset_inclusive);
        names.push_back(buff.c_str());

        buff = lfd->get_sequence(lreadidx, len, start_offset, end_offset_inclusive);
        reads.push_back(buff.c_str());

        /*! GGGG: extract kmers for this sequence but skip this sequence if the length is too short */
        if (len <= KLEN)
        {
            nskipped++;
            continue;
        }

        int nkmers = (len - KLEN + 1);
        kmersprocessed += nkmers;
        kmersthisbatch += nkmers;
        
        /* Calculate kmers */
        VectorKmer kmers = Kmer::getKmers(buff.c_str());
        ASSERT(kmers.size() == nkmers, "");
        size_t Nfound = buff.c_str().find('N');

        size_t j;
        for(j = 0; j < nkmers; ++j)
        {
            while (Nfound != string::npos && Nfound < j) 
                Nfound = buff.c_str().find('N', Nfound + 1);

            /* If there is an 'N', toss it */
            if (Nfound != string::npos && Nfound < j + KLEN)
                continue;  

            ASSERT(kmers[j] == Kmer(buff.c_str().c_str() + j), "");

            /*! GGGG: where all these variables come from? */
            size_t sending = PackEndsKmer(buff.c_str(), j, kmers[j], readIndex, j, outgoing,
                    readids, positions, extreads, pass, len, KLEN);

            if (sending > maxsending)
                maxsending = sending;
        }

        /*! GGGG: where is read id determined? */
        if (pass == 2)
        {   
            StoreReadName(readIndex, names[i], readNameMap);
        }

        /* Always start with next read index whether exiting or continuing the loop */
        readIndex++; 
        
        if (maxsending * bytesperentry >= memthreshold || (kmersthisbatch + len) * bytesperentry >= MAX_ALLTOALL_MEM)
        { 
            /* Start with next read */
            nreads = i + 1; 
            if (pass == 2)
            { 
                startReadIndex = readIndex;
            }
            break;
        }            
    } /*! GGGG: end of for all the local sequences */

  return nreads;
}

/////////////////////////////////////////////
// ProcessFiles                            //
/////////////////////////////////////////////

size_t ProcessFiles(FastaData& lfd, int pass, double& cardinality, ReadId& readIndex, std::unordered_map<ReadId, std::string>& readNameMap)
{
    /*! GGGG: include bloom filter source code */
    struct bloom * bm = NULL;
    int exchangeAndCountPass = pass;

    /* communication 
    MAX_ALLTOALL_MEM communication buffer initial size tuned for dibella 
    It's tunable if needed */
    Buffer scratch1 = initBuffer(MAX_ALLTOALL_MEM);
    Buffer scratch2 = initBuffer(MAX_ALLTOALL_MEM);

    Buffer scratch1 = initBuffer(MAX_ALLTOALL_MEM);
    Buffer scratch2 = initBuffer(MAX_ALLTOALL_MEM);
    
    /*! Initialize bloom filter */
    if(pass == 1)
    {
        /*! GGGG: random seed never used */
        bm = (struct bloom*) malloc(sizeof(struct bloom));

        const double fp_probability = 0.05;
        assert(cardinality < 1L<<32);
        bloom_init(bm, cardinality, fp_probability);

        if(myrank == 0)
        {
            std::cout << "First pass: Table size is: " << bm->bits << " bits, " << ((double)bm->bits)/8/1024/1024 << " MB" << endl;
            std::cout << "First pass: Optimal number of hash functions is : " << bm->hashes << endl;
        }

        LOGF("Initialized bloom filter with %lld bits and %d hash functions\n", (lld) bm->bits, (int) bm->hashes);
    }

    VectorVectorKmer  outgoing(nprocs);
    VectorVectorReadId readids(nprocs);
    VectorVectorChar  extreads(nprocs);
    VectorVectorPos  positions(nprocs);
        
    VectorKmer mykmers;
    VectorChar myreads;
    VectorPos  mypositions;
    VectorReadId myreadids;

    double t01 = MPI_Wtime();
    double totproctime = 0, totpack = 0, totexch = 0;
    int moreToExchange = 0;
    size_t offset = 0;
    size_t nreads = lfd->local_count();

    /* Extract kmers and counts from read sequences (seqs) */
    do {
        double texchstart = MPI_Wtime();

        /*! GGGG: Parse'n'pack, no-op if nreads == 0 */
        offset = ParseNPack(lfd, outgoing, readids, positions, readIndex, extreads, readNameMap, exchangeAndCountPass, offset);

        double tpack = MPI_Wtime() - texchstart;
        totpack += tpack;

        if (offset == nreads)
            offset = 0;

        /* Outgoing arrays will be all empty, shouldn't crush */
        double texch = ExchangePass(outgoing, readids, positions, /* extquals,*/ extreads, mykmers, myreadids, mypositions, /*myquals, myreads,*/ exchangeAndCountPass, scratch1, scratch2); 

        totexch += texch;
        DBG("Exchanged and received %lld %0.3f sec\n", (lld) mykmers.size(), texch);

        if (exchangeAndCountPass == 2)
        {
            ASSERT(mykmers.size() == myreadids.size(), "");
            ASSERT(mykmers.size() == mypositions.size(), "");
        }
        else
        {
            ASSERT(myreadids.size() == 0, "");
            ASSERT(mypositions.size() == 0, "");
        }

        /* we might still receive data even if we didn't send any */
        DealWithInMemoryData(mykmers, exchangeAndCountPass, bm, myreadids, mypositions);

        /*! GGGG: when this is the case? */
        moreToExchange = offset < nreads;

        mykmers.clear();
        myreads.clear();
        mypositions.clear();
        myreadids.clear();

        double proctime = MPI_Wtime() - texchstart - tpack - texch;
        totproctime += proctime;

        DBG("Processed (%lld): remainingToExchange = %lld %0.3f sec\n", (lld) mykmers.size(), (lld) nreads - offset, proctime);
        DBG("Checking global state: morereads = %d moreToExchange = %d moreFiles = %d\n", /* morereads, */ moreToExchange /*, moreFiles */);

        // CHECK_MPI(MPI_Allreduce(moreflags, allmore2go, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
        // DBG("Got global state: allmorereads = %d allmoreToExchange = %d allmoreFiles = %d\n", allmorereads, allmoreToExchange, allmoreFiles);

        double now = MPI_Wtime();

        if (myrank == 0 && !(exchanges % 30))
        {
            cout << __FUNCTION__ << " pass "     << pass << ": "
                //  << " active ranks morereads: "  << allmorereads
                //  << " moreToExchange: "          << allmoreToExchange
                //  << " moreFiles: "               << allmoreFiles  
                //  << ", rank "                    << myrank 
                //  << " morereads: "               << morereads 
                 << " moreToExchange: "          << moreToExchange
            cout << " tpackime: "                << std::fixed << std::setprecision( 3 ) << tpack
                 << " exchange_time: "           << std::fixed << std::setprecision( 3 ) << texch
                 << " proctimeime: "             << std::fixed << std::setprecision( 3 ) << proctime
                 << " elapsed: "                 << std::fixed << std::setprecision( 3 ) << now - t01
                 << endl;
        }

        LOGF("Exchange timings pack: %0.3f exch: %0.3f process: %0.3f elapsed: %0.3f\n", tpack, texch, proctime, now - t01);
    } while (moreToExchange);

    double t02 = MPI_Wtime();

    double tots[4], gtots[4] = {0.0, 0.0, 0.0, 0.0};

    tots[0] = totpack;
    tots[1] = totexch;
    tots[2] = totproctime;
    tots[3] = t02 - t01;

    LOGF("Process Total times: pack: %0.3f exch: %0.3f process: %0.3f\n", /* pfqTime, */ totpack, totexch, totproctime);

    CHECK_MPI(MPI_Reduce(&tots, &gtots, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));

    if (myrank == 0)
    {
        int nranks;
        CHECK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &nranks));
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for packing reads is "    << (gtots[0] / nranks) << ", myelapsed " << tots[0] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for exchanging reads is " << (gtots[1] / nranks) << ", myelapsed " << tots[1] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for processing reads is " << (gtots[2] / nranks) << ", myelapsed " << tots[2] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for elapsed is "          << (gtots[3] / nranks) << ", myelapsed " << tots[3] << endl;
    }
    
    if(myrank == 0)
    {
        cout << __FUNCTION__ << " pass " << pass << ": Read/distributed/processed reads in " << t02 - t01 << " seconds" << endl;
    }

    LOGF("Finished pass %d, freeing bloom and other memory - kmercounts: %lld entries\n", pass, (lld) kmercounts->size());

    if (bm)
    {
        LOGF("Freeing Bloom filter\n");
        bloom_free(bm);
        free(bm);
        bm = NULL;
    }

    freeBuffer(scratch1);
    freeBuffer(scratch2);

    if (exchangeAndCountPass == 2)
        countTotalKmersAndCleanHash();

    return nreads;
}

/////////////////////////////////////////////
// InsertIntoHLL                           //
/////////////////////////////////////////////

typedef struct
{
    double duration, parsingTime, getKmerTime, lexKmerTime, hllTime, hhTime;
    double last;
} MoreHLLTimers;

inline double getDuration(MoreHLLTimers &t)
{
    double delta = t.last;
    t.last = MPI_Wtime();
    return t.last - delta;
}

MoreHLLTimers InsertIntoHLL(string& myread, HyperLogLog& hll, uint64_t& found, bool extraTimers = false)
{
    MoreHLLTimers t;

    if (extraTimers)
    {
      memset(&t, 0, sizeof(MoreHLLTimers));
      getDuration(t);
      t.duration = t.last;
    }

    /* Otherwise size_t being unsigned will underflow */
    if(found >= KLEN) 
    {
        if (extraTimers) t.parsingTime += getDuration(t);

        /* Calculate all the kmers */
        std::vector<Kmer> kmers = Kmer::getKmers(myread); 
        if (extraTimers) t.getKmerTime += getDuration(t);

        ASSERT(kmers.size() >= found-KLEN+1, "");
        size_t Nfound = myread.find('N');
        
        for(size_t j = 0; j < found-KLEN+1; ++j)
        {
            ASSERT(kmers[j] == Kmer(myread.c_str() + j), "");
            while (Nfound!=std::string::npos && Nfound < j) Nfound = myread.find('N', Nfound+1);

            if (Nfound!=std::string::npos && Nfound < j+KLEN)
                continue;	/* if there is an 'N', toss it */
            Kmer &mykmer = kmers[j];

            if (extraTimers)
                t.parsingTime += getDuration(t);

            Kmer lexsmall =  mykmer.rep();
            if (extraTimers)
                t.lexKmerTime += getDuration(t);

            hll.add((const char*) lexsmall.getBytes(), lexsmall.getNumBytes());
            if (extraTimers)
                t.hllTime += getDuration(t);
        }
    }

    if (extraTimers)
    {
      t.last = t.duration;
      t.duration = getDuration(t);
    }

    return t;
}

/////////////////////////////////////////////
// ProudlyParallelCardinalityEstimate      //
/////////////////////////////////////////////

void ProudlyParallelCardinalityEstimate(FastaData& lfd, double& cardinality)
{
	HyperLogLog hll(12);
    int64_t numreads = lfd->local_count();

    /*! GGGG: so where do I get these? */
    uint64_t start_offset, end_offset_inclusive;
    uint64_t len;

    for (uint64_t lreadidx = 0; lreadidx < numreads; ++lreadidx)
    {
        /*! GGGG: loading sequence string in buff */
        char* myread = lfd->get_sequence(lseq_idx, len, start_offset, end_offset_inclusive);
		MoreHLLTimers mt = InsertIntoHLL(myread.c_str(), hll, len, true);
    }
	
    // LOGF("HLL timings: reads %lld, duration %0.4f, parsing %0.4f, getKmer %0.4f, lexKmer %0.4f, thll %0.4f, hhTime %0.4f\n", 
    // (lld) numreads, mt.duration, mt.parsingTime, mt.getKmerTime, mt.lexKmerTime, mt.thll, mt.hhTime);
    LOGF("My cardinality before reduction: %f\n", hll.estimate());

	/* Using MPI_UNSIGNED_CHAR because MPI_MAX is not allowed on MPI_BYTE */
	int count = hll.M.size();
	CHECK_MPI(MPI_Allreduce(MPI_IN_PLACE, hll.M.data(), count, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD));
	CHECK_MPI(MPI_Allreduce(MPI_IN_PLACE, &numreads, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD));

	cardinality = hll.estimate();

	if(myrank == 0)
    {
		cout << __FUNCTION__ << ": Embarrassingly parallel k-mer count estimate is " << cardinality << endl;
		cout << __FUNCTION__ << ": Total reads processed over all processors is " << readsprocessed << endl;
        ADD_DIAG("%f", "cardinality", cardinality);
        ADD_DIAG("%lld", "total_reads", (lld) readsprocessed);
	}
    
    SLOG("%s: total cardinality %f\n", __FUNCTION__, cardinality);

	MPI_Barrier(MPI_COMM_WORLD);

    /* Assume a balanced distribution */
	cardinality /= static_cast<double>(nprocs);	

    /* 10% benefit of doubt */
	cardinality *= 1.1;
    LOGF("Adjusted per-process cardinality: %f\n", cardinality);
}