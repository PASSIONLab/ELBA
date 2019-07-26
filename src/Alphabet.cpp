// Created by Saliya Ekanayake on 1/25/19.

#include <algorithm>
#include <cctype>
#include <cassert>
#include <iostream>
#include "../include/Alphabet.hpp"


/*!
 * Sometimes it is not possible two differentiate two closely related amino
 * acids, therefore we have the special cases:
 * asparagine/aspartic acid - asx - B
 * glutamine/glutamic acid - glx - Z
 * http://www.cryst.bbk.ac.uk/education/AminoAcid/the_twenty.html
 */
const char* Alphabet::protein = "ARNDCQEGHILKMFPSTWYVBZX*";
const char* Alphabet::dna = "ACGT";


Alphabet::Alphabet(Alphabet::type t) {
  switch (t){
    case PROTEIN:
      init(protein);
      size = 24;
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

  uchar code = 0;
  for (const char &c : alphabet){
    if(char_to_code[c] == no_code){
      char_to_code[c] = code;
      code_to_char[code] = static_cast<uchar>(c);
      ++code;
    }
  }
}



