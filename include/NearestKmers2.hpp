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
    char base;
    ushort sub_idx_in_base;
    short dist_to_root;

    MinSub(){}

    MinSub(char base, ushort sub_idx_in_base, short dist_to_base)
    : base(base), sub_idx_in_base(sub_idx_in_base), dist_to_root(dist_to_base){}

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
    uint64_t kmer_code;
    std::string kmer_str;
    std::unordered_set<ushort> free_idxs;
    Kmer *root_kmer = nullptr;
    ushort distance_to_root = 0;

    Kmer(std::string str, Alphabet &alph, pisa::ScoreMatrix &score_mat,
         std::unordered_set<ushort> free_idxs) :
        Kmer(std::move(str),alph, score_mat,std::move(free_idxs), this, 0){}

    Kmer(std::string str, Alphabet &alph, pisa::ScoreMatrix &score_mat,
         std::unordered_set<ushort> free_idxs, Kmer *root_kmer,
         ushort distance_to_root)
        : kmer_str(std::move(str)), free_idxs(std::move(free_idxs)),
          root_kmer(root_kmer), distance_to_root(distance_to_root) {

      kmer_code = 0;
      ushort base = alph.size;
      for (char cap_c : kmer_str) {
        /*! Efficient than using pow() */
        if (cap_c > 96 && cap_c < 123) {
          // small case character, so make it uppercase.
          cap_c = cap_c - 32;
        }
        kmer_code = kmer_code * base + alph.char_to_code[cap_c];
      }

      bool operator()(Kmer& k1, Kmer& k2) const {
        return k1.distance_to_root < k2.distance_to_root;
      }
    }

    ~Kmer(){}
  };

  class NearestKmers2 {
  public:
    NearestKmers2(Alphabet& alph, pisa::ScoreMatrix& sm);
    void print_sorted_sm();

//    std::vector<Kmer> find_nearest_kmers(Kmer root, ushort m);

    std::vector<pisa::Kmer> find_sub_kmers(Kmer root, ushort m);

  private:
    void populate_sorted_sm(Alphabet& alph, ScoreMatrix& sm);
    void explore(Kmer& p, minmax::MinMaxHeap<Kmer>& mmheap, Kmer& root, ushort m);
    void create_new_sub_kmer(Kmer& p, std::priority_queue<MinSub>& minheap,
        minmax::MinMaxHeap<Kmer>& mmheap, bool pop_max,
        ushort free_idx, ushort sub_idx, Kmer& root);
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
