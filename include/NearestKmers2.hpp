// Created by Saliya Ekanayake on 2019-10-03.

#ifndef LBL_PISA_NEARESTKMERS2_HPP
#define LBL_PISA_NEARESTKMERS2_HPP

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "Alphabet.hpp"
#include "ScoreMat.hpp"
#include "MinMaxHeap.hpp"

namespace pisa {
  typedef std::vector<std::vector<std::pair<short, char>>> SortedSM_T;

  /*!
   * Captures a triplet in the min-heap denoting a single index substitution.
   * A @b MinSub happens on a k-mer and the @dist_to_root captures the
   * accumulated substitution cost from a root k-mer to the corresponding k-mer
   * that this @b MinSub is performed plus the cost of this substitution giving
   * the total substitution cost to root k-mer if this @b MinSub were to perform.
   */
  struct MinSub{
    ushort free_idx;
    char base_at_free_idx;
    ushort sub_idx_in_base;
    short dist_to_root;

    MinSub(){}

    MinSub(ushort free_idx, char base_at_free_idx, ushort sub_idx_in_base, short dist_to_base)
    : free_idx(free_idx), base_at_free_idx(base_at_free_idx),
      sub_idx_in_base(sub_idx_in_base), dist_to_root(dist_to_base){}

    /*!
     * A greater than comparator for two @b MinSub instances.
     * This is useful in creating a min-heap out of @b MinSub instances.
     * @param ms1
     * @param ms2
     * @return
     */
    bool operator()(MinSub& ms1, MinSub& ms2) const {
      return ms1.dist_to_root > ms2.dist_to_root;
    }

  };


  struct Kmer {
  private:
    uint64_t kmer_code;
    std::string kmer_str;
    std::unordered_set<ushort> free_idxs;
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
    Kmer(std::string str, Alphabet &alph)
        : kmer_str(std::move(str)) {

      ushort free_idx = 0;
      size_t len = kmer_str.length();
      for (size_t i = 0; i < len; ++i){
        free_idxs.insert(free_idx++);
      }
      update_kmer_code(alph);
    }

    short dist_to_root() const {
      return dist2r;
    }

    Kmer substitute(ushort free_idx, char substitute_base, short dist_to_root, Alphabet& alph){
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

    char& operator[](size_t idx){
      return kmer_str[idx];
    }

    const std::unordered_set<ushort>& get_free_idxs(){
      return free_idxs;
    }


    ~Kmer(){}
  };

  class NearestKmers2 {
  public:
    NearestKmers2(Alphabet& alph, pisa::ScoreMatrix& sm);
    void print_sorted_sm();
    std::vector<pisa::Kmer> find_sub_kmers(pisa::Kmer& root, ushort m);

  private:
    void populate_sorted_sm(Alphabet& alph, ScoreMatrix& sm);
    void explore(Kmer& p,
                 minmax::MinMaxHeap<pisa::Kmer, std::vector<pisa::Kmer>,
                     pisa::Kmer>& mmheap, Kmer& root, ushort m);
    void create_new_sub_kmer(Kmer& p, std::priority_queue<pisa::MinSub, std::vector<pisa::MinSub>, pisa::MinSub>& minheap,
                             minmax::MinMaxHeap<pisa::Kmer, std::vector<pisa::Kmer>, pisa::Kmer>& mmheap, bool pop_max,
                             MinSub& ms);
    // Note, this includes the diagonal.
    SortedSM_T sorted_sm;
    Alphabet alph;
    pisa::ScoreMatrix sm;
  };


/*  static char const BLOSUM62_DATA[24*24] = {
    //A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
      4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4,
      -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4,
      -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4,
      -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4,
      0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4,
      -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4,
      -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
      0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4,
      -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4,
      -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4,
      -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4,
      -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4,
      -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4,
      -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4,
      -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4,
      1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4,
      0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4,
      -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4,
      -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4,
      0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4,
      -2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4,
      -1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4,
      0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4,
      -4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1
  };

  static const char* BLOSUM62_ALPH = "ARNDCQEGHILKMFPSTWYVBZX*";*/
}

#endif //LBL_PISA_NEARESTKMERS2_HPP
