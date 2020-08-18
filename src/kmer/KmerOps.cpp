//
// Created by Saliya Ekanayake on 2019-10-15
// Modified by Giulia Guidi on 2020-07-23
//

#include <numeric>

#include "../../include/kmer/KmerOps.hpp"
#include "../../include/HyperLogLog.hpp"

extern "C" {
#ifdef HIPMER_BLOOM64
#include "../libbloom/bloom64.h"
#else
#include "../libbloom/bloom.h"
#endif
}


#define MEGA 1000000.0
#define MILLION 1000000
#define COUNT_THRESHOLD 300000
#define COUNT_THRESHOLD_HIGH 30300000
#define HIGH_BIN 100
#define HIGH_NUM_BINS ((COUNT_THRESHOLD_HIGH-COUNT_THRESHOLD)/HIGH_BIN)

int ERR_THRESHOLD;

/*! GGGG: define this type somewhere
 *  local_kmers is not used anymore in dibella */
/*! GGGG: this contains the read id and pos information in some fency array
 * Data as I need them might aready be in line 1647 */
KmerCountsType *kmercounts = NULL;

int nprocs;
int myrank;
int64_t nonerrorkmers;
int64_t kmersprocessed;
int64_t readsprocessed;

int minKmerFreq = 2;
int maxKmerFreq = MAX_NUM_READS;

using namespace std;

namespace dibella
{

/////////////////////////////////////////////
// KmerInfo Class                          //
///////////////////////////////////////////// 

READIDS newReadIdList() {
	READIDS toreturn = *(new READIDS);
	ASSERT(nullReadId == 0, "Read ID lists are being initialized to 0, despite the reserved null read ID being non-zero"); // TODO
	std::memset(&toreturn, 0, MAX_NUM_READS*sizeof(ReadId));
	return toreturn;
}

POSITIONS newPositionsList() {
	POSITIONS toreturn = *(new POSITIONS);
	ASSERT(initPos == 0, "Position lists are being initialized to 0, despite the reserved null position being non-zero..."); // TODO
	std::memset(&toreturn, 0, MAX_NUM_READS*sizeof(PosInRead));
	return toreturn;
}

class KmerInfo {
public:
    //typedef array<char,2> TwoChar;
    //typedef unsigned int ReadId;
private:
    Kmer kmer;
    //TwoChar quals,seqs;
    ReadId readId;
    PosInRead position;
public:
    KmerInfo() {}
    KmerInfo(Kmer k): kmer(k), readId((ReadId) nullReadId), position( (PosInRead) initPos ) {}
    KmerInfo(Kmer k, ReadId r, PosInRead p): kmer(k), readId(r), position(p) { }
    //KmerInfo(Kmer k, ReadId r): kmer(k), quals(), seqs(), readId(r) {}
    //KmerInfo(Kmer k, TwoChar q, TwoChar s, ReadId r): kmer(k), quals(q), seqs(s), readId(r) {}
    KmerInfo(const KmerInfo &copy) {
        kmer = copy.kmer;
        //quals = copy.quals;
        //seqs = copy.seqs;
        readId = copy.readId;
        position = copy.position;
    }
    const Kmer& getKmer() const {
        return kmer;
    }

//     int write(GZIP_FILE f) {
//         int count = GZIP_FWRITE(this, sizeof(*this), 1, f);
// #ifndef NO_GZIP
//         if (count != sizeof(*this)*1) { DIE("There was a problem writing the kmerInfo file! %s\n", strerror(errno)); }
// #else
//         if (count != 1) { DIE("There was a problem writing the kmerInfo file! %s\n", strerror(errno)); }
// #endif
//         return count;
//     }
//     int read(GZIP_FILE f) {
//         int count = GZIP_FREAD(this, sizeof(*this), 1, f);
// #ifndef NO_GZIP
//         if (count != sizeof(*this)*1 && !GZIP_EOF(f)) { DIE("There was a problem reading the kmerInfo file! %s\n", strerror(errno)); }
// #else
//         if (count != 1 && ! feof(f)) { DIE("There was a problem reading the kmerInfo file! %s\n", strerror(errno)); }
// #endif
//         return count;
//     }

    /* Returns true if in bloom, does not modify */
    bool checkBloom(struct bloom *bm)
    {
        MPI_Pcontrol(1,"BloomFilter");
        bool inBloom = bloom_check(bm, kmer.getBytes(), kmer.getNumBytes()) == 1;
        MPI_Pcontrol(-1,"BloomFilter");
        return inBloom;	
    }

    /* Returns true if in bloom, inserts if not */
    bool checkBloomAndRemember(struct bloom *bm)
    {
        bool inBloom = checkBloom(bm);
        if (!inBloom) {
            MPI_Pcontrol(1,"BloomFilter");
            bloom_add(bm, kmer.getBytes(), kmer.getNumBytes());
            MPI_Pcontrol(-1,"BloomFilter");
        }
        return inBloom;
    }
    // returns true when kmer is already in bloom, if in bloom, inserts into map, if andCount, increments map count
    bool checkBloomAndInsert(struct bloom *bm, bool andCount) {
        bool inBloom = checkBloomAndRemember(bm);

        if (inBloom) {
            MPI_Pcontrol(1,"InsertOrUpdate");
            auto got = kmercounts->find(kmer.getArray());  // kmercounts is a global variable
            if(got == kmercounts->end())
            {
#ifdef KHASH
                kmercounts->insert(kmer.getArray(), make_tuple(newReadIdList(), newPositionsList(), 0));
#else
                kmercounts->insert(make_pair(kmer.getArray(), make_tuple(newReadIdList(), newPositionsList(), 0)));
#endif
                if (andCount) includeCount(false);
            } else {
            		if (andCount) { includeCount(got); }
            }
            MPI_Pcontrol(-1,"InsertOrUpdate");
        }
        return inBloom;
    }

    void updateReadIds(KmerCountsType::iterator got) {
#ifdef KHASH
        READIDS reads = get<0>(*got);  // ::value returns (const valtype_t &) but ::* returns (valtype_t &), which can be changed
        POSITIONS& positions = get<1>(*got);
#else
        READIDS& reads = get<0>(got->second);
        POSITIONS& positions = get<1>(got->second);
#endif
        ASSERT(readId > nullReadId,"");

        // never add duplicates, also currently doesn't support more than 1 positions per read ID
        int index;
        for (index = 0; index < maxKmerFreq && reads[index] > nullReadId; index++) {
			if (reads[index] == readId) return;
		}
        // if the loop finishes without returning, the index is set to the next open space or there are no open spaces
        if (index >= maxKmerFreq || reads[index] > nullReadId) return;

        ASSERT(reads[index] == nullReadId, "reads[index] does not equal expected value of nullReadId");
        reads[index] = readId;

        positions[index] = position;
    }

    bool includeCount(bool doStoreReadId) {
        auto got = kmercounts->find(kmer.getArray());  // kmercounts is a global variable
        if ( doStoreReadId && (got != kmercounts->end()) ) { updateReadIds(got); }
        return includeCount(got);
    }

    bool includeCount(KmerCountsType::iterator got) {
        MPI_Pcontrol(1,"HashTable");
        bool inserted = false;
        if(got != kmercounts->end()) // don't count anything else
        {
			// count the kmer in mercount
#ifdef KHASH
			++(get<2>(*got));  // ::value returns (const valtype_t &) but ::* returns (valtype_t &), which can be changed
#else
			++(get<2>(got->second)); // increment the counter regardless of quality extensions
#endif
            inserted = true;
        }
        MPI_Pcontrol(-1,"HashTable");
        return inserted;
    }
};

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

    // LOGF("my hashsize = %lld, nonerrorkmers = %lld, maxcount = %lld\n", (lld) hashsize, (lld) nonerrorkmers, (lld) maxcount);

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
            KmerInfo ki(mykmers[i]);
            ki.checkBloomAndInsert(bm, false);
        }
    }
    else
    {
        ASSERT(pass == 2, "");
        size_t count = mykmers.size();
        for(size_t i = 0; i < count; ++i)
        {
            // std::cout << "myreadids[i] in DealWithIt " << myreadids[i] << std::endl;
            KmerInfo ki(mykmers[i], myreadids[i], mypositions[i]);
			ASSERT(!bm, "");
			ki.includeCount(true);
        }
    }
}

/////////////////////////////////////////////
// ExchangePass                            //
/////////////////////////////////////////////

double ExchangePass(VectorVectorKmer& outgoing, VectorVectorReadId& readids, VectorVectorPos& positions, VectorVectorChar& extreads,
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

    // DBG("totsend = %lld totrecv = %lld\n", (lld) totsend, (lld) totrecv);

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
        // extquals[i].clear();
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

    // serial_printf("KmerMatch:%s exchange iteration %d pass %d: sent min %lld bytes, sent max %lld bytes, recv min %lld bytes, recv max %lld bytes, in min %.3f s, max %.3f s\n",
    //     __FUNCTION__, iter, pass, global_mins[SND], global_maxs[SND], global_mins[RCV], global_maxs[RCV], global_min_time, global_max_time);

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

    // DBG("DeleteAll: recvcount = %lld, sendct = %lld\n", (lld) recvcnt, (lld) sendcnt);
    DeleteAll(rdispls, sdispls, recvcnt, sendcnt);

    // iter++;
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
        int pass, int lastCountedBase)
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
size_t ParseNPack(FastaData* lfd, VectorVectorKmer& outgoing, VectorVectorReadId& readids, VectorVectorPos& positions, 
    ReadId& startReadIndex, VectorVectorChar& extreads, std::unordered_map<ReadId, std::string>& readNameMap, int pass, size_t offset)
{
    uint64_t start_offset, end_offset_inclusive;
    ushort len;

    /* offset left over from orginal piece of codes */
    size_t nreads = lfd->local_count();
    size_t maxsending = 0, kmersthisbatch = 0, nskipped = 0;
    size_t bytesperkmer  = Kmer::numBytes();
    size_t bytesperentry = bytesperkmer + 4;
    size_t memthreshold  = (MAX_ALLTOALL_MEM/nprocs) * 2;

    /*! GGGG: where this startReadIndex come from? */
    ReadId readIndex = startReadIndex;

    vector<string> names;
    vector<string> reads;

    char* buff;

    /*! GGGG: make sure name are consistent with ids */
    for (uint64_t lreadidx = offset; lreadidx < nreads; ++lreadidx)
    {
        buff = lfd->get_sequence_id(lreadidx, len, start_offset, end_offset_inclusive);
        string myname = string(buff);
        names.push_back(myname);

        buff = lfd->get_sequence(lreadidx, len, start_offset, end_offset_inclusive);
        string myseq = string(buff);
        reads.push_back(myseq);


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
        VectorKmer kmers = Kmer::getKmers(myseq);
        ASSERT(kmers.size() == nkmers, "");
        size_t Nfound = myseq.find('N');

        size_t j;
        for(j = 0; j < nkmers; ++j)
        {
            while (Nfound != string::npos && Nfound < j) 
                Nfound = myseq.find('N', Nfound + 1);

            /* If there is an 'N', toss it */
            if (Nfound != string::npos && Nfound < j + KLEN)
                continue;  

            ASSERT(kmers[j] == Kmer(myseq.c_str() + j), "");

            /*! GGGG: where all these variables come from? */
            size_t sending = PackEndsKmer(myseq, j, kmers[j], readIndex, j, outgoing,
                    readids, positions, extreads, pass, len);

            if (sending > maxsending)
                maxsending = sending;
        }

        /*! GGGG: where is read id determined? */
        if (pass == 2)
        {   
            StoreReadName(readIndex, names[lreadidx], readNameMap);
        }

        /* Always start with next read index whether exiting or continuing the loop */
        readIndex++; 
        
        if (maxsending * bytesperentry >= memthreshold || (kmersthisbatch + len) * bytesperentry >= MAX_ALLTOALL_MEM)
        { 
            /* Start with next read */
            nreads = lreadidx + 1; 
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

size_t ProcessFiles(FastaData* lfd, int pass, double& cardinality, ReadId& readIndex, std::unordered_map<ReadId, std::string>& readNameMap)
{
    /*! GGGG: include bloom filter source code */
    struct bloom * bm = NULL;
    int exchangeAndCountPass = pass;

    /* communication 
    MAX_ALLTOALL_MEM communication buffer initial size tuned for dibella 
    It's tunable if needed */
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
            std::cout << __FUNCTION__ << ": First pass: Table size is: " << bm->bits << " bits, " << ((double)bm->bits)/8/1024/1024 << " MB" << endl;
            std::cout << __FUNCTION__ << ": First pass: Optimal number of hash functions is : " << bm->hashes << endl;
        }
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

    int exchanges = 0;

    /* Extract kmers and counts from read sequences (seqs) */
    do {
        double texchstart = MPI_Wtime();

        /*! GGGG: Parse'n'pack, no-op if nreads == 0 */
        offset = ParseNPack(lfd, outgoing, readids, positions, readIndex, extreads, readNameMap, exchangeAndCountPass, offset);

        double tpack = MPI_Wtime() - texchstart;
        totpack += tpack;

        /* GGGG: change this logic --too confusing; this temporary change might introduce a bug on big data set */
        // if (offset == nreads)
        //     offset = 0;

        /* Outgoing arrays will be all empty, shouldn't crush */
        double texch = ExchangePass(outgoing, readids, positions, /* extquals,*/ extreads, mykmers, myreadids, mypositions, /*myquals, myreads,*/ exchangeAndCountPass, scratch1, scratch2); 

        totexch += texch;
        // DBG("Exchanged and received %lld %0.3f sec\n", (lld) mykmers.size(), texch);

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

        std::cout << __FUNCTION__ << ": offset " << offset << std::endl;
        std::cout << __FUNCTION__ << ": nreads " << nreads << std::endl;

        moreToExchange = offset < nreads;

        mykmers.clear();
        myreads.clear();
        mypositions.clear();
        myreadids.clear();

        double proctime = MPI_Wtime() - texchstart - tpack - texch;
        totproctime += proctime;

        // DBG("Processed (%lld): remainingToExchange = %lld %0.3f sec\n", (lld) mykmers.size(), (lld) nreads - offset, proctime);
        // DBG("Checking global state: morereads = %d moreToExchange = %d moreFiles = %d\n", /* morereads, */ moreToExchange /*, moreFiles */);

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
                 << " moreToExchange: "          << moreToExchange;
            cout << " tpackime: "                << std::fixed << std::setprecision(3) << tpack
                 << " exchange_time: "           << std::fixed << std::setprecision(3) << texch
                 << " proctimeime: "             << std::fixed << std::setprecision(3) << totproctime
                 << " elapsed: "                 << std::fixed << std::setprecision(3) << now - t01
                 << endl;
        }

        // LOGF("Exchange timings pack: %0.3f exch: %0.3f process: %0.3f elapsed: %0.3f\n", tpack, texch, proctime, now - t01);
    } while (moreToExchange);

    double t02 = MPI_Wtime();

    double tots[4], gtots[4] = {0.0, 0.0, 0.0, 0.0};

    tots[0] = totpack;
    tots[1] = totexch;
    tots[2] = totproctime;
    tots[3] = t02 - t01;

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

    if (bm)
    {
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

void ProudlyParallelCardinalityEstimate(FastaData* lfd, double& cardinality)
{
	HyperLogLog hll(12);
    int64_t numreads = lfd->local_count();

    /*! GGGG: so where do I get these? */
    uint64_t start_offset, end_offset_inclusive, found;
    ushort len;

    for (uint64_t lreadidx = 0; lreadidx < numreads; ++lreadidx)
    {
        /*! GGGG: loading sequence string in buff */
        char* myread = lfd->get_sequence(lreadidx, len, start_offset, end_offset_inclusive);

        std::string mystr = std::string(myread);
        found = len;
        
		MoreHLLTimers mt = InsertIntoHLL(mystr, hll, found, true);
    }
    
    // LOGF("HLL timings: reads %lld, duration %0.4f, parsing %0.4f, getKmer %0.4f, lexKmer %0.4f, thll %0.4f, hhTime %0.4f\n", 
    // (lld) numreads, mt.duration, mt.parsingTime, mt.getKmerTime, mt.lexKmerTime, mt.thll, mt.hhTime);
    // LOGF("My cardinality before reduction: %f\n", hll.estimate());

	/* Using MPI_UNSIGNED_CHAR because MPI_MAX is not allowed on MPI_BYTE */
	int count = hll.M.size();
	CHECK_MPI(MPI_Allreduce(MPI_IN_PLACE, hll.M.data(), count, MPI_UNSIGNED_CHAR, MPI_MAX, MPI_COMM_WORLD));
	CHECK_MPI(MPI_Allreduce(MPI_IN_PLACE, &numreads, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD));

	cardinality = hll.estimate();

	if(myrank == 0)
    {
		cout << __FUNCTION__ << ": Embarrassingly parallel k-mer count estimate is " << cardinality << endl;
		cout << __FUNCTION__ << ": Total reads processed over all processors is " << readsprocessed << endl;
	}
    
    // SLOG("%s: total cardinality %f\n", __FUNCTION__, cardinality);

	MPI_Barrier(MPI_COMM_WORLD);

    /* Assume a balanced distribution */
	cardinality /= static_cast<double>(nprocs);	

    /* 10% benefit of doubt */
	cardinality *= 1.1;

	if(myrank == 0)
    {
		cout << __FUNCTION__ << ": Adjusted per-process cardinality: " << cardinality << endl;
	}
}

PSpMat<POSITIONS>::MPI_DCCols KmerOps::generate_A(uint64_t seq_count,
      std::shared_ptr<DistributedFastaData> &dfd, ushort k, ushort s,
      Alphabet &alph, const std::shared_ptr<ParallelOps> &parops,
      const std::shared_ptr<TimePod> &tp) /*, std::unordered_set<Kmer, Kmer>& local_kmers) */
  {

  char *buff;
  ushort len;
  uint64_t start_offset, end_offset_inclusive;

  /* typedef std::vector<uint64_t> uvec_64; */
  uvec_64 lrow_ids, lcol_ids;
  std::vector<PosInRead> lvals;

  uint64_t offset = dfd->global_start_idx();
  FastaData *lfd  = dfd->lfd();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: Cardinality estimation
/////////////////////////////////////////////////////////////////////////////////////////////////////////

  tp->times["start_kmerop:gen_A:loop_add_kmers()"] = std::chrono::system_clock::now();

  /*! GGGG: cardinality estimate */
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  unsigned int readsxproc = 0;
  int myrank = parops->world_proc_rank;

  double cardinality;
  double tstart = MPI_Wtime(); 

  Kmer::set_k(KLEN);
  
  /* Doesn't update kmersprocessed yet (but updates readsprocessed */
  ProudlyParallelCardinalityEstimate(lfd, cardinality); 
	readsxproc = readsprocessed / nprocs;

  double tcardinalitye = MPI_Wtime() - tstart;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: First pass k-mer counter
/////////////////////////////////////////////////////////////////////////////////////////////////////////

  /*! Prepare for first pass over input data with k-mer extraction */
  ASSERT(kmercounts == NULL, "");

  /*! GGGG: what is this? */
  kmercounts = new KmerCountsType();

  /*! Assume we have at least 10x depth of kmers to store and a minimum 64k */
  uint64_t reserve = cardinality * 0.1 + 64000; 
  /*! This is the maximum supported by khash */
  if (reserve >= 4294967291ul) reserve = 4294967290ul; // 

  /*! GGGG: define this macros */
  // LOGF("Reserving %lld entries in VectorMap for cardinality %lld\n", (lld) reserve, (lld) cardinality);
  kmercounts->reserve(reserve);
  // DBG("Reserved kmercounts\n");

  /* Initialize readNameMap for storing ReadID -> names/tags of reads */
  /* GGGG: define ReadId type */
  std::unordered_map<ReadId, std::string>* readNameMap = new std::unordered_map<ReadId, std::string>();

  /*! GGGG: I don't what the original one, I wnt the new one with consecutive entries; also I only need the first one; it's gonna be incremented later in ParseNPack */
  uint64_t GlobalReadOffset = dfd->g_seq_offsets[parops->world_proc_rank];  
  ReadId myReadStartIndex = GlobalReadOffset;

  /*! GGGG: let's extract the function (I'll separate later once I understood what's going on) */
  /*! GGGG: functions in KmerCounter.cpp */
  /*  Determine final hash-table entries using bloom filter */
  int nreads = ProcessFiles(lfd, 1, cardinality, myReadStartIndex, *readNameMap);//, readids);

  double firstpasstime = MPI_Wtime() - tstart;

  /*! GGGG: TODO print using PASTIS way */
  if (myrank == 0)
  {
    std::cout << "Reads: " << nreads << std::endl;
    std::cout << "First input data pass, elapsed time: " << firstpasstime << std::endl;
  }   

  tstart = MPI_Wtime();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: Second pass k-mer counter        
/////////////////////////////////////////////////////////////////////////////////////////////////////////

  /*! GGGG: prepare for second pass over the input data, extracting k-mers with read ID's, names, etc. */
  /* Exchange number of reads-per-processor to calculate read indices */
  uint64_t sndReadCounts[nprocs];

  for (int i = 0; i < nprocs; i++)
  {
    sndReadCounts[i] = nreads;
  }

  uint64_t recvReadCounts[nprocs];

  CHECK_MPI(MPI_Alltoall(sndReadCounts, 1, MPI_UINT64_T, recvReadCounts, 1, MPI_UINT64_T, MPI_COMM_WORLD));
  
  myReadStartIndex = 1;
  for (int i = 0; i < myrank; i++)
  {
    myReadStartIndex += recvReadCounts[i];
  }

#ifdef DEBUG
    if(myrank == nprocs - 1)
      for (int i = 0; i < nprocs; i++)
        LOGF("recvReadCounts[%lld] = %lld\n", recvReadCounts[i]);
#endif

  CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

  /* Initialize end read ranges */
  ReadId readRanges[nprocs];
  ReadId last = 0;

  for (int i = 0; i < nprocs; i++)
  {
    readRanges[i] = last + recvReadCounts[i];
    last += recvReadCounts[i];
  }

  // DBG("My read range is [%lld - %lld]\n", (myrank==0? 1 : readRanges[myrank-1]+1), readRanges[myrank]);

  /* Second pass */
  ProcessFiles(lfd, 2, cardinality, myReadStartIndex, *readNameMap);//, readids);

  double timesecondpass = MPI_Wtime() - tstart;
  // serial_printf("%s: 2nd input data pass, elapsed time: %0.3f s\n", __FUNCTION__, timesecondpass);
  tstart = MPI_Wtime();

  int64_t sendbuf = kmercounts->size(); 
  int64_t recvbuf, totcount, maxcount;

  CHECK_MPI(MPI_Exscan(&sendbuf, &recvbuf, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
  CHECK_MPI(MPI_Allreduce(&sendbuf, &totcount, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
  CHECK_MPI(MPI_Allreduce(&sendbuf, &maxcount, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD));

  int64_t totkmersprocessed;

  CHECK_MPI(MPI_Reduce(&kmersprocessed, &totkmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));

  double timeloadimbalance = MPI_Wtime() - tstart;

  if(myrank  == 0)
  {
    cout << __FUNCTION__ << ": Total number of stored k-mers: " << totcount << endl;

    double imbalance = static_cast<double>(nprocs * maxcount) / static_cast<double>(totcount);  
    
    cout << __FUNCTION__ << ": Load imbalance for final k-mer counts: " << imbalance << endl;
    cout << __FUNCTION__ << ": CardinalityEstimate " << static_cast<double>(totkmersprocessed) / (MEGA * max((tcardinalitye), 0.001) * nprocs) << " MEGA k-mers per sec/proc in " << (tcardinalitye) << " seconds"<< endl;
    cout << __FUNCTION__ << ": Bloom filter + hash table (key) initialization " << static_cast<double>(totkmersprocessed) / (MEGA * max((firstpasstime),0.001) * nprocs) << " MEGA k-mers per sec/proc in " << (firstpasstime) << " seconds" << endl;
    cout << __FUNCTION__ << ": Hash table (value) initialization  " << static_cast<double>(totkmersprocessed) / (MEGA * max((timesecondpass),0.001) * nprocs) << " MEGA k-mers per sec/proc in " << (timesecondpass) << " seconds" << endl;
  }

  // serial_printf("%s: Total time computing load imbalance: %0.3f s\n", __FUNCTION__, timeloadimbalance);
  CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
 
  tstart = MPI_Wtime();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: Build matrix A
/////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*! GGGG: Once k-mers are consolidated in a single global location, 
     *  the other processors don’t need to know the ids of kmers they sent off to other processors. */

    std::unordered_map<Kmer::MERARR, uint64_t>* kmerIdMap = new std::unordered_map<Kmer::MERARR, uint64_t>();

    uint64_t kmerid = 0;
    for(auto itr = kmercounts->begin(); itr != kmercounts->end(); ++itr)
    {
        
        /*! GGGG: TODO assing ids to local kmers here, they need to be consecutive on procs */
        /*! kmer string */
        Kmer::MERARR key = itr->first;
        Kmer mykmer(key);

        /*! GGGG: run tests, read idx should be consistent now */
        READIDS  readids = get<0>(itr->second);
        POSITIONS values = get<1>(itr->second);

        for(int j = 0; j < readids.size(); j++)
        {
            if(readids[j] != 0)
            {
            /*!  GGGG: I need kmer id here 
             *   This temp solution only works on 1 node */
                lcol_ids.push_back(kmerid);
                lrow_ids.push_back(readids[j]);
                lvals.push_back(values[j]);

                std::cout << "k " << kmerid     << " " << mykmer << std::endl;
                std::cout << "r " << readids[j] << std::endl;
                std::cout << "v " << values[j]  << std::endl;
            }
        }
        kmerid++;
    }

    assert(kmerid == kmercounts->size());
    exit(0); 

#ifndef NDEBUG
  {
    std::string title = "Local matrix info:";
    std::string msg   = "lrow_ids row_size: " + 
                     std::to_string(lrow_ids.size())
                     + " lcol_ids row_size: " +
                     std::to_string(lcol_ids.size())
                     + " lvals row_size: " + 
                     std::to_string(lvals.size());
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

  assert(lrow_ids.size() == lcol_ids.size() && lcol_ids.size() == lvals.size());

//   /*! Create distributed sparse matrix of sequence x kmers */
//   FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, parops->grid);
//   FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, parops->grid);
//   FullyDistVec<uint64_t, POSITIONS> dvals(lvals, parops->grid);

//   uint64_t nrows = seq_count;
//   /*! Columns of the matrix are direct maps to kmers identified
//    * by their |alphabet| base number. E.g. for proteins this is
//    * base 20, so the direct map has to be 20^k in size. */

//   /*! GGGG: this is the reliable k-mer space */
//   auto ncols = static_cast<uint64_t>(pow(alph.size, k));
//   tp->times["start_kmerop:gen_A:spMatA()"] = std::chrono::system_clock::now();
//   PSpMat<POSITIONS>::MPI_DCCols A(nrows, ncols, drows, dcols, dvals, false);
//   tp->times["end_kmerop:gen_A:spMatA()"]   = std::chrono::system_clock::now();
  
//   return A;
  }
}
