// Created by Saliya Ekanayake on 2019-10-03.

#include <iostream>
#include <cassert>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include "../include/NearestKmers2.hpp"

void populate_blosum62(pisa::SortedSM_T& sorted_sm){
  for (size_t i = 0; i < 24; ++i){
    size_t row_offset = i*24;
    char row_c = pisa::BLOSUM62_ALPH[i];
    short self_score = pisa::BLOSUM62_DATA[row_offset];
    for (size_t j = 1; j < 24; ++j){
      char col_c = pisa::BLOSUM62_ALPH[j];
      short score = pisa::BLOSUM62_DATA[row_offset+j];

    }
  }

}

pisa::NearestKmers2::NearestKmers2()  {
  populate_blosum62(sorted_sm);
}

int main(int argc, char** argv){
  boost::uuids::random_generator gen;
  boost::uuids::uuid id = gen();

  std::cout << id << '\n';
  return 0;
}