#ifndef DIBELLA_KMER_HPP
#define DIBELLA_KMER_HPP

#include <stdio.h>
#include <stdint.h>
#include <cassert>
#include <cstring>
#include <string>
#include <array>
#include <vector>
#include <functional>
#include <cstdint>

/*! GGGG: add these files from diBELLA */
#include "../Common.h"
#include "../HashFuncs.h"

/*! GGGG: Kmer.cpp/hpp come from diBELLA */
/* Short description: 
 *  - Store kmer strings by using 2 bits per base instead of 8 
 *  - Easily return reverse complements of kmers, e.g. TTGG -> CCAA
 *  - Easily compare kmers
 *  - Provide hash of kmers
 *  - Get last and next kmer, e.g. ACGT -> CGTT or ACGT -> AACGT
 *  */
#define N_LONGS (MAX_KMER_SIZE/32)
#define N_BYTES (MAX_KMER_SIZE/4)

class Kmer {
  public:
    typedef std::array<uint64_t, N_LONGS> MERARR;
    typedef std::array<uint8_t, N_BYTES> BYTEARR;
    
    Kmer();
    Kmer(const Kmer& o);
    explicit Kmer(const char *s);

    /* This is like a shadow constructor (to avoid accidental signature match with the existing constructor) */
    void copyDataFrom(uint8_t *  mybytes)	
    {
      memcpy(longs.data(), mybytes, sizeof(uint64_t) * (N_LONGS));
    }

    explicit Kmer(const MERARR & arr)
    {
    	std::memcpy (longs.data(), arr.data(), sizeof(uint64_t) * (N_LONGS));
    }
    static std::vector<Kmer> getKmers(std::string seq);
    
    Kmer& operator=(const Kmer& o);  
    void set_deleted();
    bool operator<(const Kmer& o) const;
    bool operator==(const Kmer& o) const;
    
    bool operator!=(const Kmer& o) const
    {
      return !(*this == o);
    }
    
    void set_kmer(const char *s);
    uint64_t hash() const;
    
    Kmer twin() const;
    /* ABAB: return the smaller of itself (lexicographically) or its reversed-complement (i.e. twin) */
    Kmer rep() const; 
    Kmer getLink(const size_t index) const;
    Kmer forwardBase(const char b) const;
    Kmer backwardBase(const char b) const;
    std::string getBinary() const;  
    void toString(char * s) const;
    std::string toString() const;
    
    void copyDataInto(void * pointer) const
    {
    	  memcpy(pointer, longs.data(), sizeof(uint64_t) * (N_LONGS));
    }
    
    /* ABAB: return the raw data packed in an std::array
     * this preserves the lexicographical order on k-mers
     * i.e. A.toString() < B.toString <=> A.getArray() < B.getArray()
     */
    const MERARR &getArray() const {
          return longs;
    }
    const uint8_t *getBytes() const {
          return bytes.data();
    }
    int getNumBytes() const {
          return N_BYTES;
    }
    
    /* Returns true for completely identical k-mers as well as k-mers that only differ at the last base */
    bool equalUpToLastBase(const Kmer & rhs);	

    static void set_k(unsigned int _k);
    static constexpr size_t numBytes() { 
    	  return (sizeof(uint64_t) * (N_LONGS));
    }
    
    static const unsigned int MAX_K = MAX_KMER_SIZE;
    static unsigned int k;
    
  private:
    static unsigned int k_bytes;
    static unsigned int k_modmask;
    
    /* Data fields */
    union {
      MERARR  longs;
      BYTEARR bytes;
    };
    
    // Unions are very useful for low-level programming tasks that involve writing to the same memory area 
    // but at different portions of the allocated memory space, for instance:
    //		union item {
    //			// The item is 16-bits
    //			short theItem;
    //			// In little-endian lo accesses the low 8-bits -
    //			// hi, the upper 8-bits
    //			struct { char lo; char hi; } portions;
    //		};
    //  item tItem;
    //  tItem.theItem = 0xBEAD;
    //  tItem.portions.lo = 0xEF; // The item now equals 0xBEEF

};

struct KmerHash {
  size_t operator()(const Kmer &km) const {
    return km.hash();
  }
};

/* Specialization of std::Hash */
namespace std
{
  template<> struct hash<Kmer>
  {
    typedef std::size_t result_type;
    result_type operator()(Kmer const& km) const
    {
      return km.hash();
    }
  }
  
  template<> struct hash<Kmer::MERARR>
  {
    typedef std::size_t result_type;
    result_type operator()(const Kmer::MERARR & km) const
    {
      return MurmurHash3_x64_64((const void*)km.data(),sizeof(Kmer::MERARR));
    }
  };
};

inline std::ostream& operator<<(std::ostream& out, const Kmer& k){
    return out << k.toString();
};

/*! GGGG: MAX_NUM_READS defined in CMakeFiles.txt */
/*! Currently records 1 position per (k-mer, read) pair */
typedef std::array<PosInRead, MAX_NUM_READS> POSITIONS;
typedef std::array<ReadId,    MAX_NUM_READS> READIDS;

typedef tuple<READIDS, POSITIONS, int> KmerCountType;
typedef pair<Kmer::MERARR, KmerCountType>  KmerValue;

/*! GGGG: might need to modify this */
/*! GGGG: import vector map */
typedef VectorMap<Kmer::MERARR, KmerCountType, std::hash<Kmer::MERARR>, std::less<Kmer::MERARR>, std::equal_to<Kmer::MERARR>> KmerCountsType;

#endif // DIBELLA_KMER_HPP
