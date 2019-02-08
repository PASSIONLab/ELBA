// Created by Saliya Ekanayake on 1/22/19.


#ifndef LBL_DAL_ALPHABET_HPP
#define LBL_DAL_ALPHABET_HPP

#include <string>
#include "Types.hpp"

/*!
 * Data structure to hold and manipulate info about the sequence alphabet.
 */
struct Alphabet {
  static const char* protein;
  static const char* dna;

  static const ushort capacity = 128;

  void init(const std::string &letters);

  enum type {PROTEIN, DNA};
  uchar char_to_code[capacity];
  uchar code_to_char[capacity];
  ushort size;

  explicit Alphabet(type t);
};

#endif //LBL_DAL_ALPHABET_HPP
