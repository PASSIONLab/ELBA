// Created by Saliya Ekanayake on 1/22/19.


#ifndef LBL_DAL_ALPHABET_HPP
#define LBL_DAL_ALPHABET_HPP

#include <string>

/*!
 * Data structure to hold and manipulate info about the sequence alphabet.
 */
struct Alphabet {
  static const char* protein;
  static const char* dna;

  static const unsigned capacity = 128;

  void init(const std::string &letters);

  enum type {PROTEIN, DNA};
  unsigned char char_to_code[10];
  unsigned char code_to_char[10];

  Alphabet(type t);
};

#endif //LBL_DAL_ALPHABET_HPP
