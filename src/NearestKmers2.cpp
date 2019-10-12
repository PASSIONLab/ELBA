// Created by Saliya Ekanayake on 2019-10-03.

#include <iostream>
#include <cassert>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include "../include/NearestKmers2.hpp"

/*!
 * Sorts the given scoring matrix by substitution cost
 * to create @b sorted_sm, which is a 2D vector.
 * The first dimension is any character @b c in the alphabet.
 * The second dimension is sorted such that for any two
 * indices @b i and @b j in [0, |@b alphabet|),
 * sorted_sm[c][i] less than or equal to sorted_sm[c][j]
 * if (sm.score(c, c) - sm.score(c, alph[i])) less than or equal to (sm.score(c, c) - sm.score(c, alph[j]))

 *
 * @param alph
 * @param sm
 */
void pisa::NearestKmers2::populate_sorted_sm(Alphabet& alph, pisa::ScoreMatrix& sm) {
  for (size_t i = 0; i < alph.size; ++i){
    char row_c = alph.letters[i];
    short self_score = sm.score(row_c, row_c);
    for (size_t j = 0; j < alph.size; ++j){
      char col_c = alph.letters[j];
      short score = sm.score(row_c, col_c);
      sorted_sm[row_c].emplace_back(std::make_pair((self_score - score), col_c)); // TODO - fix emplace back
    }
    std::sort(sorted_sm[row_c].begin(), sorted_sm[row_c].end(),
              [](auto& p1, auto& p2){
                return p1.first < p2.first;
              });
  }
}

std::vector<pisa::Kmer> pisa::NearestKmers2::find_sub_kmers(pisa::Kmer root, ushort m){

}




/*std::vector<pisa::Kmer> pisa::NearestKmers2::find_nearest_kmers(pisa::Kmer root,
                                                                ushort m) {

  // Let's generate the first level of nearest kmers
  ushort free_idx_count = root.free_idxs.size();
  std::vector<size_t> current_idxs(free_idx_count, 1);
  std::priority_queue<pisa::MinSub, std::vector<pisa::MinSub>, pisa::MinSub> min_pq;

  // Put the first entries (except diagonal or self) to heap
//  for (ushort free_idx : root.free_idxs){
//    MinSub ms(sorted_sm[])
//  }

  std::vector<pisa::Kmer> ret;
  return ret;

}*/

void pisa::NearestKmers2::print_sorted_sm() {
  for (int i = 0; i < 24; ++i){
    char row_c = alph.letters[i];
    for (auto p : sorted_sm[row_c]){
      std::cout << "(" << p.first << "," << p.second << ") ";
    }
    std::cout << std::endl;
  }
}

pisa::NearestKmers2::NearestKmers2(Alphabet& alph, pisa::ScoreMatrix& sm)
: sorted_sm(alph.max_char+1), sm(sm), alph(alph)  {
  populate_sorted_sm(alph, sm);
}

void
pisa::NearestKmers2::explore(
    pisa::Kmer& p,
    minmax::MinMaxHeap<pisa::Kmer>& mmheap,
    pisa::Kmer& root, ushort m) {

}

void
pisa::NearestKmers2::create_new_sub_kmer(
    Kmer& p, std::priority_queue<MinSub>& minheap,
    minmax::MinMaxHeap<Kmer>& mmheap, bool pop_max,
    ushort free_idx, ushort sub_idx, Kmer& root)
{
  minheap.pop();
//  ushort dist2root = sorted_sm[]

}

int main(int argc, char** argv){
//  boost::uuids::random_generator gen;
//  boost::uuids::uuid id = gen();
//
//  std::cout << id << '\n';
  /*std::vector<std::vector<int>> test(10);
  std::cout << (int)'A' << std::endl;
  std::cout << (int)'Z' << std::endl;
  std::cout << (int)'*' << std::endl;
  std::cout << test[5].size() << std::endl;*/

  Alphabet alph(Alphabet::PROTEIN);
  pisa::Blosum62 sm;

  pisa::NearestKmers2 nk2(alph, sm);
  nk2.print_sorted_sm();

  int M = 5; // 5 Nearest k-mers

  std::unordered_set<ushort> free_idxs;
  std::string kmer_str{"ACT"};
  for (ushort i = 0; i < kmer_str.length(); ++i) {
    free_idxs.insert(i);
  }
  // Root k-mer
  pisa::Kmer root(kmer_str, alph, sm, free_idxs);

  std::vector<pisa::Kmer> nearest_kmers = nk2.find_nearest_kmers(root, M);


  std::priority_queue<pisa::MinSub, std::vector<pisa::MinSub>, pisa::MinSub> pq;
  pisa::MinSub ms1('A', 2, 10);
  pisa::MinSub ms2('B', 2, 5);
  pq.push(ms1);
  pq.push(ms2);

  while(!pq.empty()){
    std::cout << pq.top().dist_to_root << std::endl;
    pq.pop();
  }

  return 0;
}