
// Created by Saliya Ekanayake on 2019-09-12.

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <stack>
#include "../include/Alphabet.hpp"
#include "../include/ScoreMat.hpp"


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
  }

  ~Kmer(){}
};

struct KmerComp {
  bool operator()(Kmer *&lhs, Kmer *&rhs) const {
    return lhs->distance_to_root < rhs->distance_to_root;
  }
};

std::vector<Kmer> findNNearestKmers(int num_kmers, Kmer &root,
                                    Alphabet &alph,
                                    pisa::ScoreMatrix &score_mat) {
  std::vector<Kmer> nearest_kmers;
  std::unordered_set<uint64_t> explored_kmers;
  // k-mers that get thrown away from the heap
  // We'll need to this to check if a k-mer we pick from the
  // to_explore_kmers is to be discarded or not.
  std::unordered_set<uint64_t> expunged_kmers;
  std::priority_queue<Kmer *, std::vector<Kmer *>, KmerComp> heap;
  std::stack<Kmer *> to_explore_kmers;
  to_explore_kmers.push(&root);


  /* TODO (Saliya): unique_ptr instead of raw pointers would be nice */
  while (!to_explore_kmers.empty()) {
    int to_explore_size = to_explore_kmers.size();
    for (int explore_idx = 0; explore_idx < to_explore_size; ++explore_idx) {
      Kmer *current = to_explore_kmers.top();
      to_explore_kmers.pop();
      if (expunged_kmers.find(current->kmer_code) != expunged_kmers.end()) {
        // This has already been removed from the heap, so no need to explore.
        continue;
      }

//      if (explored_kmers.find(current) != explored_kmers.end()){
//        // This has already been explored in a previous round, so skip it.
//        continue;
//      }
      for (ushort free_idx : current->free_idxs) {
        std::unordered_set<ushort> child_free_idxs = current->free_idxs;
        child_free_idxs.erase(free_idx);
        for (const char &c : alph.letters) {
          if (current->kmer_str[free_idx] == c) {
            continue;
          }
          // Swap char in free_idx of current with c to create a new child current.
          std::string child_kmer_str = current->kmer_str;
          child_kmer_str[free_idx] = c;

          short score_root_current
              = score_mat.score(current->root_kmer->kmer_str[free_idx],
                                current->kmer_str[free_idx]);
          short score_current_child
              = score_mat.score(current->kmer_str[free_idx], c);
          ushort dist_to_root =
              std::abs(score_root_current - score_current_child)
              + current->distance_to_root;

          // Check if this child kmer survives the heap.
          // The conditions are,
          // 1. if heap size < num_kmers, child can go in (i.e. skips the "if")
          // 2. if heap full (size == num_kmers) then
          //    2.1 if child's distance to root is  > heap's top's dist to root
          //        then child cannot go in (i.e. goes into "if")
          //    OR
          //    2.2 if the child and heap's top both have the same dist to root
          //        but heap's top has more free idxs (i.e. less mutated)
          //        then child can't go in (i.e. goes into "if")
          if (!heap.empty()){
            Kmer *top = heap.top();
            if (heap.size() == num_kmers &&
                (top->distance_to_root < dist_to_root ||
                 (top->distance_to_root == dist_to_root &&
                  top->free_idxs.size() > child_free_idxs.size()))) {
              continue;
            }
            // put the current top element to expunged list
            expunged_kmers.insert(top->kmer_code);
            heap.pop();
          }


          // This child needs to go into the heap
          Kmer *child = new Kmer(child_kmer_str, alph, score_mat,
                                 child_free_idxs,
                                 current->root_kmer, dist_to_root);


//          delete (top); // delete all later
          heap.push(child);
        }
      }
      explored_kmers.insert(current->kmer_code);
    } // END of all kmers in the to_explore_list for this phase

    // At this point we should sync the to_explore_kmers and heap
    // such that they to_explore_kmers contains a subset of the heap and
    // they are sorted by the distance to root in ascending order.
    // Note. it's possible for our heap to contain kmers that were in the
    // just explored set of kmers, so when pushing into the to_explore_kmer
    // stack we'll check against the explored_kmer set so we'll not push
    // already explored kmers back into the to_explore_stack.
    // This way we can run our outerloop until no more kmers to explore.
    auto heap_copy = heap;
    while(!heap_copy.empty()){
      Kmer* top = heap_copy.top();
      heap_copy.pop();
      if (explored_kmers.find(top->kmer_code) != explored_kmers.end()){
        continue;
      }
      // This will ensure the sorted order.
      to_explore_kmers.push(top);
    }
  }

  for (int i = 0; i < num_kmers; ++i){
    Kmer* top = heap.top();
    nearest_kmers.push_back(*top);
    if (top->root_kmer != top){
      // TODO bug here
//      delete(top);
    }
  }

  // Clean explored kmers
  // TODO fix here

  std::cout << root.kmer_code << std::endl;


  return nearest_kmers;
}

int main(int argc, char **argv) {
  int N = 5; // 5 Nearest k-mers
  Alphabet alph(Alphabet::PROTEIN);
  pisa::Blosum62 bsum62;

  std::unordered_set<ushort> free_idxs;
  std::string kmer_str{"ACT"};
  for (ushort i = 0; i < kmer_str.length(); ++i) {
    free_idxs.insert(i);
  }
  // Root k-mer
  Kmer kmer(kmer_str, alph, bsum62, free_idxs);
  std::vector<Kmer> nearest_kmers = findNNearestKmers(N, kmer, alph, bsum62);
  std::cout << ">>Nearest kmers\n";
  for (auto& kmer : nearest_kmers){
    std::cout << kmer.kmer_str << std::endl;
  }
//  std::string test("ABC");
//  const char* alph = "DDC";
//  std::string alph_str(alph);
//  test[0] = alph_str[0];
//  std::cout << test << std::endl;
//  std::cout << alph_str << std::endl;
}