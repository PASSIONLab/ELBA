#ifndef _DIBELLA_DEFINES_H_
#define _DIBELLA_DEFINES_H_

#include <stdlib.h> 
#include <array>
#include <map>
#include <vector>
#include <utility>
#include <tuple>
#include <iostream>
#include <string>

#ifndef MAX_ALLTOALL_MEM
#define MAX_ALLTOALL_MEM (128*1024*1024)  /* 128 MB */
//#define MAX_ALLTOALL_MEM (12*1024*1024) /*  12 MB */
#endif

typedef long long int lld;
typedef unsigned long long int llu;

#ifdef DEBUG
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::terminate(); \
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#ifndef __FILENAME__
#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#endif

#define EXIT_FUNC(x) do { common_exit(x); } while (0)
#ifndef  __UPC_VERSION__
#define MYTHREAD get_rank()
#define THREADS get_num_ranks()
#endif

#ifndef DEFAULT_READ_LEN
  #define DEFAULT_READ_LEN 70000
#endif

/* Maximum number of alignments is read_length - kmer_length + 1 */
#ifndef MAX_ALIGN
#define MAX_ALIGN DEFAULT_READ_LEN - MAX_KMER_SIZE + 1
#endif

#ifndef ONE_KB
  #define ONE_KB 1024L
  #define ONE_MB 1048576L
  #define ONE_GB 1073741824L
#endif

/* for printing contigs to file */
#ifndef SEGMENT_LENGTH
  #define SEGMENT_LENGTH 51
#endif

#ifndef LIB_NAME_LEN // TODO check if this value is appropriate for long reads (was directly copied from hipmer for loadfq)
#define LIB_NAME_LEN 15
#endif

#ifdef __linux__ // TODO temporary - this is not the best way to do this
#define VMFILE_PATH "/dev/shm/"
#else // TODO
#define VMFILE_PATH ""
#endif

#define INIT_READIDS_FNAME "InitReadIdx.txt"
#define OVERLAP_OUT_FNAME  "MyOverlap.out"
#define ALGNMNT_OUT_FNAME  "MyAlignment.out"

/*
 * A fixed-sized vector of read ID's,
 * identified by integers.
 */
typedef uint64_t ReadId;
const ReadId nullReadId = 0;

#ifdef TIGHT_READS
typedef uint16_t PosInRead; // read lengths up to 65,535 bps
#else
typedef uint32_t PosInRead; // read lengths > 65,535 bps
#endif
const PosInRead initPos = 0; // assumes positions will only be processed if ReadId is non-null

//typedef std::pair< std::pair<ReadId, ReadId>, std::pair<PosInRead, PosInRead> > ReadOverlapPair;
typedef std::pair<PosInRead, PosInRead> PosPair;
typedef std::map< std::pair<ReadId, ReadId>, std::vector<PosPair> > ReadOverlapMap;

struct ReadOverlapPair {
	ReadId readId1, readId2;
	PosInRead posInRead1, posInRead2;
	ReadOverlapPair(): readId1(nullReadId), readId2(nullReadId), posInRead1(0), posInRead2(0) {}
	ReadOverlapPair(ReadId r, ReadId s, PosInRead p, PosInRead q): readId1(r), readId2(s), posInRead1(p), posInRead2(q) {}
} ;

inline std::ostream & operator<<(std::ostream & str, const ReadOverlapPair& pair) {
	return (str << "[" << pair.readId1 << ", " << pair.readId2 << ": " << pair.posInRead1 << ", " << pair.posInRead2 <<"]");
}

#define GET_KMER_PACKED_LEN(k) ((k + 3) / 4)

#endif /* _DIBELLA_DEFINES_H_ */