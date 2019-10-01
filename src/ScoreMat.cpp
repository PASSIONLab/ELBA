// Created by Saliya Ekanayake on 2019-07-23.

#include "../include/ScoreMat.hpp"
#include "../include/Alphabet.hpp"


pisa::ScoreMatrix::ScoreMatrix(ushort alph_size) :
  alph_size(alph_size){
}

pisa::Blosum62::Blosum62() : ScoreMatrix(24) {
  const char *alph = Alphabet::protein;
  for (int i = 0; i < alph_size; ++i){

    char ci = alph[i];
    ushort score_offset = ci * row_size;
    ushort data_offset = i * alph_size;

    std::vector<char> *subs = new std::vector<char>(alph_size);
    short max_score = SHRT_MIN;
    short max_score_count = 0;

    for (int j = 0; j < alph_size; ++j){
      char cj = alph[j];
      char s = data[data_offset + j];

      _score[score_offset + cj] = s;

      if (ci != cj){
        if (s > max_score){
          max_score_count = 0;
          max_score = s;
          (*subs)[max_score_count] = cj;
          ++max_score_count;
        } else if (s == max_score){
          (*subs)[max_score_count] = cj;
          ++max_score_count;
        }
      }
    }
    (*subs)[max_score_count] = ci;
    ++max_score_count;
    subs->erase(subs->begin()+max_score_count, subs->end());

#ifndef NDEBUG
    std::ostringstream vts;
    // Convert all but the last element to avoid a trailing ","
    std::copy(subs->begin(), subs->end()-1,
              std::ostream_iterator<char>(vts, ", "));

    // Now add the last element with no delimiter
    vts << subs->back();
    std::cout << "subs for " << ci << " " << subs->size() << " " << vts.str() << std::endl;
#endif

    base_to_subtitutes[ci] = subs;
  }

#ifndef NDEBUG
  /*! For debug purposes */
  for (int i = 0; i < alph_size; ++i){
    char ci = alph[i];
    ushort offset = ci * row_size;
    for (int j = 0; j < alph_size; ++j){
      char cj = alph[j];
      char s = _score[offset + cj];
      assert(s == data[i * alph_size + j]);
      std::cout << (short)s << " ";
    }
    std::cout << std::endl;
  }
#endif



}

pisa::Blosum62::~Blosum62(){
  for (std::pair<char, std::vector<char>*> p : base_to_subtitutes){
    free(p.second);
  }
}
