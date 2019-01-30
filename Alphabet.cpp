// Created by Saliya Ekanayake on 1/25/19.

#include <algorithm>
#include <cctype>
#include <cassert>
#include <iostream>
#include "Alphabet.hpp"

const char* Alphabet::protein = "ACDEFGHIKLMNPQRSTVWY";
const char* Alphabet::dna = "ACGT";


Alphabet::Alphabet(Alphabet::type t) {
  switch (t){
    case PROTEIN:
      init(protein);
      size = 20;
      break;
    case DNA:
      init(dna);
      size = 4;
      break;
  }
}

void Alphabet::init(const std::string &alphabet) {
  unsigned no_code = capacity;
  std::fill_n(char_to_code, capacity, no_code);

  unsigned code = 0;
  for (const char &c : alphabet){
    if(char_to_code[c] == no_code){
      char_to_code[c] = static_cast<unsigned char>(code);
      code_to_char[code] = (unsigned char) c;
      ++code;
    }
  }
}



