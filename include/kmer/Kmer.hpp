// Created by Saliya Ekanayake on 10/17/19.

#ifndef DIBELLA_KMER_HPP
#define DIBELLA_KMER_HPP
#include <set>
#include <vector>
#include "../Alphabet.hpp"
//#include "ScoreMat.hpp"
//#include "MinMaxHeap.hpp"

namespace dibella
{
  struct Kmer
  {
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
}

/*! GGGG: MAX_NUM_READS defined in CMakeFiles.txt */
/*! Currently records 1 position per (k-mer, read) pair */
typedef std::array<PosInRead, MAX_NUM_READS> POSITIONS;
typedef std::array<ReadId, MAX_NUM_READS> READIDS;

typedef tuple<READIDS, POSITIONS, int>  KmerCountType;
typedef pair<Kmer::kmer_str, KmerCountType> KmerValue;

/*! GGGG: might need to modify this */
/*! GGGG: import vector map */
typedef VectorMap<Kmer::MERARR, KmerCountType, std::hash<Kmer::MERARR>, std::less<Kmer::MERARR>, std::equal_to<Kmer::MERARR>> KmerCountsType;

#endif //DIBELLA_KMER_HPP
