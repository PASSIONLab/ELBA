// Created by Saliya Ekanayake on 10/17/19.

#ifndef DIBELLA_KMER_HPP
#define DIBELLA_KMER_HPP
#include <set>
#include <vector>
#include "../Alphabet.hpp"
#include "../HashFuncs.h"
#include "../Defines.hpp"
#include "../VectorMap.hpp"

namespace dibella {
  struct Kmer {
  private:
    uint64_t kmer_code;
    std::string kmer_str;
    // We need to get free_idxs in reverse sorted order during k-mer generation
    // so can't use an unordered_set here.
    // Note. the order doesn't matter for the nearest k-mer generation.
    std::set<ushort, std::greater<ushort>> free_idxs;
    short dist2r = 0;

    void update_kmer_code(Alphabet& alph){
      ushort base = alph.size;
      kmer_code = 0;
      for (char cap_c : kmer_str) {
        /*! Efficient than using pow() */
        if (cap_c > 96 && cap_c < 123) {
          // small case character, so make it uppercase.
          cap_c = cap_c - 32;
        }
        kmer_code = kmer_code * base + alph.char_to_code[cap_c];
      }
    }

  public:
    Kmer(){}
    Kmer(uint64_t kmer_code, ushort k, Alphabet& alph)
    : kmer_code(kmer_code) {
      uint64_t q, r;
      ushort free_idx = 0;
      for (ushort i = 0; i < k-1; ++i) {
        q = kmer_code / alph.size;
        r = kmer_code - (q*alph.size);
        kmer_str.insert(kmer_str.begin(), alph.code_to_char[r]);
        free_idxs.insert(free_idx++);
        kmer_code = q;
      }
      kmer_str.insert(kmer_str.begin(), alph.code_to_char[kmer_code]);
      free_idxs.insert(free_idx);
    }

    Kmer(std::string str, uint64_t kmer_code, Alphabet& alph, bool fix_kcode)
    : kmer_str(std::move(str)), kmer_code(kmer_code){
      ushort free_idx = 0;
      size_t len = kmer_str.length();
      for (size_t i = 0; i < len; ++i){
        free_idxs.insert(free_idx++);
      }
      if (fix_kcode) {
        update_kmer_code(alph);
      }
    }

    Kmer(std::string str, Alphabet &alph)
        : Kmer(str, 0, alph, true) {

    }

    inline uint64_t code() const{
      return kmer_code;
    }

    inline std::string str() const{
      return kmer_str;
    }

    inline short dist_to_root() const {
      return dist2r;
    }

    Kmer substitute(ushort free_idx, char substitute_base, short dist_to_root, Alphabet& alph) const{
      Kmer subk = *this;
      subk.kmer_str[free_idx] = substitute_base;
      subk.update_kmer_code(alph);
      subk.free_idxs.erase(free_idx);
      subk.dist2r = dist_to_root;
      return subk;
    }

    bool operator()(const Kmer& k1, const Kmer& k2) const {
      return k1.dist2r < k2.dist2r;
    }

    friend std::ostream & operator << (std::ostream &out, const Kmer &k)
    {
      out << k.kmer_str << " kcode: " << k.kmer_code << " dist: " << k.dist_to_root();
      out << " free_idxs: { ";
      for (auto& free_idx : k.free_idxs){
        out << free_idx << " ";
      }
      out << "}";
      ushort hop = k.kmer_str.length() - k.free_idxs.size();
      out << " " << hop << "-hop" << std::endl;
      return out;
    }

    char operator[](size_t idx) const {
      return kmer_str[idx];
    }

    const std::set<ushort, std::greater<ushort>>& get_free_idxs() const{
      return free_idxs;
    }

    bool operator==(const Kmer& t) const
    {
      return (kmer_code == t.kmer_code);
    }

    size_t operator()(const Kmer& t) const
    {
      return t.kmer_code;
    }

    ~Kmer(){}
  };

    // GGGG: MAX_KMER_SIZE is defined in CMakeFile.txt
    /*! Currently records 1 position per (k-mer, read) pair */
    typedef std::array<PosInRead, MAX_NUM_READS> POSITIONS;
    typedef std::array<ReadId,    MAX_NUM_READS> READIDS;

    typedef tuple<READIDS, POSITIONS, int> KmerCountType;
    typedef pair<Kmer, KmerCountType>  KmerValue;

    // GGGG: might need to modify this 
    typedef VectorMap<Kmer, KmerCountType, std::hash<Kmer>, std::less<Kmer>, std::equal_to<Kmer>> KmerCountsType;
}

/////////////////////////////////////////////
// KmerInfo                                //
///////////////////////////////////////////// 

// class KmerInfo {
// public:
//     //typedef array<char,2> TwoChar;
//     //typedef unsigned int ReadId;
// private:
//     Kmer kmer;
//     //TwoChar quals,seqs;
//     ReadId readId;
//     PosInRead position;
// public:
//     KmerInfo() {}
//     KmerInfo(Kmer k): kmer(k), readId( (ReadId) nullReadId), position( (PosInRead) initPos ) {}
//     KmerInfo(Kmer k, ReadId r, PosInRead p): kmer(k), readId(r), position(p) {}
//     //KmerInfo(Kmer k, ReadId r): kmer(k), quals(), seqs(), readId(r) {}
//     //KmerInfo(Kmer k, TwoChar q, TwoChar s, ReadId r): kmer(k), quals(q), seqs(s), readId(r) {}
//     KmerInfo(const KmerInfo &copy) {
//         kmer = copy.kmer;
//         //quals = copy.quals;
//         //seqs = copy.seqs;
//         readId = copy.readId;
//         position = copy.position;
//     }
//     const Kmer& getKmer() const {
//         return kmer;
//     }
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
//     // returns true if in bloom, does not modify
//     bool checkBloom(struct bloom *bm) {
//         MPI_Pcontrol(1,"BloomFilter");
//         bool inBloom = bloom_check(bm, kmer.getBytes(), kmer.getNumBytes()) == 1;
//         MPI_Pcontrol(-1,"BloomFilter");
//         return inBloom;	
//     }
//     // returns true if in bloom, inserts if not
//     bool checkBloomAndRemember(struct bloom *bm) {
//         bool inBloom = checkBloom(bm);
//         if (!inBloom) {
//             MPI_Pcontrol(1,"BloomFilter");
//             bloom_add(bm, kmer.getBytes(), kmer.getNumBytes());
//             MPI_Pcontrol(-1,"BloomFilter");
//         }
//         return inBloom;
//     }
//     // returns true when kmer is already in bloom, if in bloom, inserts into map, if andCount, increments map count
//     bool checkBloomAndInsert(struct bloom *bm, bool andCount) {
//         bool inBloom = checkBloomAndRemember(bm);

//         if (inBloom) {
//             MPI_Pcontrol(1,"InsertOrUpdate");
//             auto got = kmercounts->find(kmer.getArray());  // kmercounts is a global variable
//             if(got == kmercounts->end())
//             {
// #ifdef KHASH
//                 kmercounts->insert(kmer.getArray(), make_tuple(newReadIdList(), newPositionsList(), 0));
// #else
//                 kmercounts->insert(make_pair(kmer.getArray(), make_tuple(newReadIdList(), newPositionsList(), 0)));
// #endif
//                 if (andCount) includeCount(false);
//             } else {
//             		if (andCount) { includeCount(got); }
//             }
//             MPI_Pcontrol(-1,"InsertOrUpdate");
//         }
//         return inBloom;
//     }

//     void updateReadIds(KmerCountsType::iterator got) {
// #ifdef KHASH
//         READIDS reads = get<0>(*got);  // ::value returns (const valtype_t &) but ::* returns (valtype_t &), which can be changed
//         POSITIONS& positions = get<1>(*got);
// #else
//         READIDS& reads = get<0>(got->second);
//         POSITIONS& positions = get<1>(got->second);
// #endif
//         ASSERT(readId > nullReadId,"");

//         // never add duplicates, also currently doesn't support more than 1 positions per read ID
//         int index;
//         for (index = 0; index < reliable_max && reads[index] > nullReadId; index++) {
// 			if (reads[index] == readId) return;
// 		}
//         // if the loop finishes without returning, the index is set to the next open space or there are no open spaces
//         if (index >= reliable_max || reads[index] > nullReadId) return;
//         ASSERT(reads[index] == nullReadId, "reads[index] does not equal expected value of nullReadId");
//         reads[index] = readId;
//         	positions[index] = position;
//     }

//     bool includeCount(bool doStoreReadId) {
//         auto got = kmercounts->find(kmer.getArray());  // kmercounts is a global variable
//         if ( doStoreReadId && (got != kmercounts->end()) ) { updateReadIds(got); }
//         return includeCount(got);
//     }

//     bool includeCount(KmerCountsType::iterator got) {
//         MPI_Pcontrol(1,"HashTable");
//         bool inserted = false;
//         if(got != kmercounts->end()) // don't count anything else
//         {
// 			// count the kmer in mercount
// #ifdef KHASH
// 			++(get<2>(*got));  // ::value returns (const valtype_t &) but ::* returns (valtype_t &), which can be changed
// #else
// 			++(get<2>(got->second)); // increment the counter regardless of quality extensions
// #endif
//             inserted = true;
//         }
//         MPI_Pcontrol(-1,"HashTable");
//         return inserted;
//     }
// };

#endif //DIBELLA_KMER_HPP
