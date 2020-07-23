/////////////////////////////////////////////
// ountTotalKmersAndCleanHash              //
///////////////////////////////////////////// 

/*! GGGG: strat from here */

/////////////////////////////////////////////
// ExchangePass                            //
/////////////////////////////////////////////

double Exchange(vector< vector< Kmer > >& outgoing, vector< vector< ReadId > >& readids, vector< vector< PosInRead > >& positions, vector< vector< array< char, 2 > > >& extseqs,
              vector< Kmer >& mykmers, vector< ReadId >& myreadids, vector< PosInRead >& mypositions, int pass, Buffer scratch1, Buffer scratch2)
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

    int * sendcnt = new int[nprocs];

    for(int i = 0; i < nprocs; ++i)
    {
        sendcnt[i] = (int) outgoing[i].size() * bytesperentry;

        if (pass == 2)
        {
            ASSERT( outgoing[i].size() == readids[i].size(), "" );
            ASSERT( outgoing[i].size() == positions[i].size(), "" );
        }
        else
        {
            ASSERT (readids[i].size() == 0, "");
            ASSERT (positions[i].size() == 0, "");
        }
    }

    int * sdispls = new int[nprocs];
    int * rdispls = new int[nprocs];
    int * recvcnt = new int[nprocs];

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
        extseqs[i].clear();
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
inline size_t FinishPackPass1(vector< vector<Kmer> > & outgoing, Kmer & kmerreal)
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
inline size_t FinishPackPass2(vector< vector<Kmer> >& outgoing, vector< vector< ReadId > >& readids, vector< vector< PosInRead > >& positions, 
    vector< vector< array< char, 2 > > >& extquals, vector< vector< array< char, 2 > > >& extseqs, Kmer& kmerreal, ReadId readId, PosInRead pos)
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

size_t PackEndsKmer(string& seq, int j, Kmer &kmerreal, ReadId readid, PosInRead pos, vector< vector<Kmer> >& outgoing,
        vector< vector<ReadId> >& readids, vector< vector<PosInRead> >& positions, vector< vector< array< char, 2 > > >& extquals,
        vector< vector< array< char, 2 > > >& extseqs, int pass, int lastCountedBase, int kmer_length)
{
    bool isCounted = lastCountedBase >= j + kmer_length;
    size_t procSendCount;

    assert(seq.size() >= j + kmer_length);
    assert(seq.substr(j, kmer_length).find('N') == std::string::npos);
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
        procSendCount = FinishPackPass2(outgoing, readids, positions, extquals, extseqs, kmerreal, readid, pos);
    }
    return procSendCount;
}

/////////////////////////////////////////////
// ProcessFiles                            //
/////////////////////////////////////////////

size_t ProcessFiles(const vector<filedata> & allfiles, int pass, double& cardinality, bool cached_io,
                        const char* base_dir, ReadId& readIndex, std::unordered_map<ReadId, std::string>& readNameMap)
{
    /*! GGGG: include bloom filter source code */
    struct bloom * bm = NULL;
    int exchangeAndCountPass = pass;

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

    vector< vector< Kmer > >  outgoing(nprocs);
    vector< vector< ReadId > > readids(nprocs);
    vector< vector< array< char, 2 > > > extquals(nprocs);
    vector< vector< array< char, 2 > > > extreads(nprocs);
    vector< vector< PosInRead> > positions(nprocs);
        
    vector< Kmer > mykmers;
    vector< ReadId > myreadids;
    vector< array< char, 2 > > myreads;
    vector< PosInRead > mypositions;

    vector< string > reads;
    vector< string > names;

    /////////////////////////////////////////////
    // ParseNPack                              //
    /////////////////////////////////////////////
    
    /*! GGGG: Parse'n'pack */
    size_t nreads = lfd->local_count();
    size_t nskipped   = 0;
    size_t maxsending = 0, kmersthisbatch = 0;
    size_t bytesperkmer  = Kmer::numBytes();
    size_t bytesperentry = bytesperkmer + 4;
    size_t memthreshold = (MAX_ALLTOALL_MEM/nprocs) * 2;

    ReadId readIndex = startReadIndex;

    for (uint64_t lseq_idx = 0; lseq_idx < lfd->local_count(); ++lseq_idx)
    {
      /*! GGGG: loading sequence id tag in buff */
      buff = lfd->get_sequence_id(lseq_idx, len, start_offset, end_offset_inclusive);
      names.push_back(buff.c_str());

      /*! GGGG: loading sequence string in buff */
      buff = lfd->get_sequence(lseq_idx, len, start_offset, end_offset_inclusive);
      reads.push_back(buff.c_str());

      /*! GGGG: calculate the max read len for each of the input files and write to marker files */)
      // writeMaxReadLengths(maxReadLengths, allfiles);

      /*! GGGG: extract kmers for this sequence */
      /* Skip this sequence if the length is too short */
      if (len <= KMER_LENGTH)
      {
          nskipped++;
          continue;
      }

      int nkmers = (len - KMER_LENGTH + 1);
      kmersprocessed += nkmers;
      kmersthisbatch += nkmers;
      
      /* Calculate kmers */
      vector< Kmer > kmers = Kmer::getKmers(buff.c_str());
      ASSERT(kmers.size() == nkmers, "");
      size_t Nfound = buff.c_str().find('N');

      size_t j;
      for(j = 0; j < nkmers; ++j)
      {
          while (Nfound != string::npos && Nfound < j) 
              Nfound = buff.c_str().find('N', Nfound + 1);

          /* If there is an 'N', toss it */
          if (Nfound != string::npos && Nfound < j + KMER_LENGTH)
              continue;  

          ASSERT(kmers[j] == Kmer(buff.c_str().c_str() + j), "");

          /*! GGGG: where all these variables come from? */
          size_t sending = PackEndsKmer(buff.c_str(), j, kmers[j], readIndex, j, outgoing,
                  readids, positions, extquals, extseqs, pass, len, KMER_LENGTH);

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
      
      if (maxsending * bytesperentry >= memoryThreshold || (kmersthisbatch + len) * bytesperentry >= MAX_ALLTOALL_MEM)
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

    /////////////////////////////////////////////
    // End of ParseNPack                       //
    /////////////////////////////////////////////

    if (offset == reads.size())
    {
        /* no need to do the swap trick as we will reuse these buffers in the next iteration */
        /*! GGGG: I might need the swap trick */ 
        reads.clear();
        offset = 0;
    }

    /* Outgoing arrays will be all empty, shouldn't crush */
    double texch = ExchangePass(outgoing, readids, positions, /* extquals,*/ extseqs, mykmers, myreadids, mypositions, /*myquals, myseqs,*/ exchangeAndCountPass, scratch1, scratch2); 

#ifdef DEBUG
    if(myrank == 0) 
        cout << "Finished Exchange pass " << exchangeAndCountPass << endl;
#endif

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
    moreToExchange = offset < seqs.size();

    mykmers.clear();
    myreadids.clear();
    mypositions.clear();
    myseqs.clear();

    double proctime = MPI_Wtime() - texchstart - tpack - texch;
    totproctime += proctime;

    DBG("Processed (%lld).  remainingToExchange = %lld %0.3f sec\n", (lld) mykmers.size(), (lld) seqs.size() - offset, proctime);
    DBG("Checking global state: moreSeqs = %d moreToExchange = %d moreFiles = %d\n", moreSeqs, moreToExchange, moreFiles);

    CHECK_MPI( MPI_Allreduce(moreflags, allmore2go, 3, MPI_INT, MPI_SUM, MPI_COMM_WORLD) );
    DBG("Got global state: allmoreSeqs = %d allmoreToExchange = %d allmoreFiles = %d\n", allmoreSeqs, allmoreToExchange, allmoreFiles);

    double now = MPI_Wtime();

    if (myrank == 0 && !(exchanges % 30))
    {
        cout << __FUNCTION__ << " pass "     << pass << ": "
             << " active ranks moreSeqs: "   << allmoreSeqs
             << " moreToExchange: "          << allmoreToExchange
             << " moreFiles: "               << allmoreFiles  
             << ", rank "                    << myrank 
             << " moreSeqs: "                << moreSeqs 
             << " moreToExchange: "          << moreToExchange
             << " moreFiles: "               << moreFiles;
        cout << " tpackime: "                << std::fixed << std::setprecision( 3 ) << tpack
             << " exchange_time: "           << std::fixed << std::setprecision( 3 ) << texch
             << " proctimeime: "             << std::fixed << std::setprecision( 3 ) << proctime
             << " elapsed: "                 << std::fixed << std::setprecision( 3 ) << now - t01
             << endl;
    }

    LOGF("Exchange timings pack: %0.3f exch: %0.3f process: %0.3f elapsed: %0.3f\n", tpack, texch, proctime, now - t01);

    t02 = MPI_Wtime();

    /*! GGGG: might not need this */
    traw = traw - pfqTime;

    double tots[6], gtots[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    tots[0] = pfqTime;
    tots[1] = tpack;
    tots[2] = texch;
    tots[3] = totproctime;
    tots[4] = traw;
    tots[5] = t02 - t01;

    LOGF("Process Total times: fastq: %0.3f pack: %0.3f exch: %0.3f process: %0.3f elapsed: %0.3f\n", pfqTime, tpack, texch, totproctime, tots[4]);

    CHECK_MPI(MPI_Reduce(&tots, &gtots, 6, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));

    if (myrank == 0)
    {
        int num_ranks;
        CHECK_MPI(MPI_Comm_size(MPI_COMM_WORLD, &num_ranks));
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for FASTQ reads is " << (gtots[0] / num_ranks)    << ", myelapsed " << tots[0] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for packing reads is " << (gtots[1] / num_ranks)  << ", myelapsed " << tots[1] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for exchanging reads is " << (gtots[2] / num_ranks) << ", myelapsed " << tots[2] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for processing reads is " << (gtots[3] / num_ranks) << ", myelapsed " << tots[3] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for other FASTQ processing is " << (gtots[4] / num_ranks) << ", myelapsed " << tots[4] << endl;
        cout << __FUNCTION__ << " pass " << pass << ": Average time taken for elapsed is " << (gtots[5] / num_ranks) << ", myelapsed " << tots[5] << endl;
    }
    
    if(myrank == 0)
    {
        cout << __FUNCTION__ << " pass " << pass << ": Read/distributed/processed reads of " << (files_itr == allfiles.end() ? " ALL files " : files_itr->filename) << " in " << t02 - t01 << " seconds" << endl;
    }

    LOGF("Finished pass %d, Freeing bloom and other memory. kmercounts: %lld entries\n", pass, (lld) kmercounts->size());

    t02 = MPI_Wtime();

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

    /////////////////////////////////////////////
    // Summary 7/22/2020                       //
    ///////////////////////////////////////////// 

    /*! GGGG: I extracted kmers from my local reads and exchange them with other processes; 
     *  I haven't extracted positions or read id yet --- I just now the tag id of my reads. 
     *  TODO: countTotalKmersAndCleanHash();
     */

    return nreads;
}
