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

  enum type {PROTEIN, DNA};
  explicit Alphabet(type t);
  void init(const std::string &letters);

  uchar char_to_code[capacity];
  uchar code_to_char[capacity];
  ushort size;
  /*! The character with the maximum code.
   * For example, in the protein alphabet this
   * would be Z as its character code is 90.
   */
  ushort max_char;
  std::string letters;
};



#endif //LBL_DAL_ALPHABET_HPP
