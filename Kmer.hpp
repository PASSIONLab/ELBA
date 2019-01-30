// Created by Saliya Ekanayake on 1/29/19.

#ifndef LBL_DAL_KMER_HPP
#define LBL_DAL_KMER_HPP

#include <cstdint>
#include <vector>
#include <cmath>
#include "Alphabet.hpp"

typedef std::vector<int64_t> vec64_t;

int count_kmers(const char* seq, int len, int start_offset,
                 int enf_offset_inclusive, int k, int s, Alphabet& alp,
                 vec64_t& lcol_ids, vec64_t& lvals){
  auto num_kmers = static_cast<int>(floor((len - k)*1.0 / s)) + 1;
  // TODO: Saliya - this can be improved using bit operators
  int64_t base = alp.size;
  int64_t kcode = 0;
  int count = 0;
  for (int i = start_offset; i <= enf_offset_inclusive; i+=s){
    kcode = 0;
    if (count == num_kmers) break;
    for (int j = i; j < i+k; ++j){
      /*! Efficient than using pow() */
      kcode = kcode * base + alp.char_to_code[*(seq + j)];
    }
    ++count;
    lcol_ids.push_back(kcode);
    lvals.push_back(i);
  }

  return num_kmers;
}



#endif //LBL_DAL_KMER_HPP
