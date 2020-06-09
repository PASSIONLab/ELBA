// Created by Saliya Ekanayake on 2019-10-03.

#include <iostream>
#include <cassert>
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
void dibella::NearestKmers2::populate_sorted_sm(Alphabet& alph, dibella::ScoreMatrix& sm) {
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

std::vector<dibella::Kmer>
    dibella::NearestKmers2::find_sub_kmers(const dibella::Kmer& root, ushort m)
{
  std::vector<dibella::Kmer> nbrs;
  minmax::MinMaxHeap<dibella::Kmer, std::vector<dibella::Kmer>, dibella::Kmer> mmheap;
  explore(root, mmheap, root, m);
  while(nbrs.size() < m){
    dibella::Kmer min_kmer = mmheap.findMin();
    nbrs.push_back(min_kmer);
    if (min_kmer.get_free_idxs().size() > 0) {
      // Explore further only if this min_kmer has 
      // free_idxs left to modify. Else, just continue
      // picking the next min_kmer from the mmheap.
      explore(min_kmer, mmheap, root, m);
    }
    mmheap.popMin();
  }
  return nbrs;
}

void dibella::NearestKmers2::print_sorted_sm() {
  for (int i = 0; i < 25; ++i){
    char row_c = alph.letters[i];
    for (auto p : sorted_sm[row_c]){
      std::cout << "(" << p.first << "," << p.second << ") ";
    }
    std::cout << std::endl;
  }
}

dibella::NearestKmers2::NearestKmers2(Alphabet& alph, dibella::ScoreMatrix& sm)
: sorted_sm(alph.max_char+1), sm(sm), alph(alph)  {
  populate_sorted_sm(alph, sm);
}

void
dibella::NearestKmers2::explore(
    const dibella::Kmer& p,
    minmax::MinMaxHeap<dibella::Kmer, std::vector<dibella::Kmer>, dibella::Kmer>& mmheap,
    const dibella::Kmer& root, ushort m) {
  std::priority_queue<dibella::MinSub, std::vector<dibella::MinSub>, dibella::MinSub> minheap;
  for (auto& free_idx : p.get_free_idxs()){
    char base_at_free_idx = p[free_idx];
    minheap.emplace(free_idx, base_at_free_idx, 1,
        p.dist_to_root()+sorted_sm[base_at_free_idx][1].first);
  }

  if (mmheap.size() < m){
    do{
      MinSub ms = minheap.top();
      create_new_sub_kmer(p, minheap, mmheap, false, ms);
    } while (mmheap.size() != m);
  } else {
    bool eliminated = false;
    do{
      MinSub ms = minheap.top();
      if(ms.dist_to_root < mmheap.findMax().dist_to_root()){
        create_new_sub_kmer(p, minheap, mmheap, true, ms);
      } else {
        eliminated = true;
      }
    } while(!eliminated);
  }
}

void
dibella::NearestKmers2::create_new_sub_kmer(
    const Kmer& p, std::priority_queue<dibella::MinSub, std::vector<dibella::MinSub>, dibella::MinSub>& minheap,
    minmax::MinMaxHeap<Kmer, std::vector<dibella::Kmer>, dibella::Kmer>& mmheap, bool pop_max, MinSub& ms)
{
  minheap.pop();
  std::pair<short, char> cost_and_subc = sorted_sm[ms.base_at_free_idx][ms.sub_idx_in_base];
  Kmer sub_kmer = p.substitute(ms.free_idx, cost_and_subc.second, ms.dist_to_root, alph);
  if (pop_max){
    mmheap.popMax();
  }
  mmheap.push(sub_kmer);
  if (ms.sub_idx_in_base+1 < alph.size){
    auto dist_to_root = (short) (p.dist_to_root() +
                            sorted_sm[ms.base_at_free_idx][ms.sub_idx_in_base+1].first);
    minheap.emplace(ms.free_idx, ms.base_at_free_idx,
        ms.sub_idx_in_base+1, dist_to_root);
  }
}

int main2(int argc, char** argv)
{
//  boost::uuids::random_generator gen;
//  boost::uuids::uuid id = gen();
//
//  std::cout << id << '\n';
  /*std::vector<std::vector<int>> test(10);
  std::cout << (int)'A' << std::endl;
  std::cout << (int)'Z' << std::endl;
  std::cout << (int)'*' << std::endl;
  std::cout << test[5].size() << std::endl;*/

  Alphabet alph(Alphabet::DNA);
  dibella::Blosum62 sm;

  dibella::NearestKmers2 nk2(alph, sm);
  nk2.print_sorted_sm();

  int M = 5; // 5 Nearest k-mers

//  std::string kmer_str{"TACTBZDP"};
  std::string kmer_str{"****"};
  // Root k-mer
  dibella::Kmer root(kmer_str, alph);
  dibella::Kmer subk = root.substitute(0, 'T', sm.dist(root[0], 'T'), alph);
  std::cout << root << std::endl;
  std::cout << subk << std::endl;
  subk = subk.substitute(2, 'G', sm.dist(subk[2], 'G'), alph);
  std::cout << subk << std::endl;

//  std::vector<dibella::Kmer> nearest_kmers = nk2.find_nearest_kmers(root, M);


  std::priority_queue<dibella::MinSub, std::vector<dibella::MinSub>, dibella::MinSub> pq;
  dibella::MinSub ms1(0, 'A', 2, 10);
  dibella::MinSub ms2(1, 'B', 2, 5);
  pq.push(ms1);
  pq.push(ms2);

  while(!pq.empty()){
    std::cout << pq.top().dist_to_root << std::endl;
    pq.pop();
  }

  std::vector<dibella::Kmer> nbrs = nk2.find_sub_kmers(root, 50);
  for (auto& kmer : nbrs){
    std::cout << kmer;
  }
  return 0;
}
