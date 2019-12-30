// Created by Saliya Ekanayake on 10/15/19.

#ifndef LBL_PISA_COMMONKMERS_HPP
#define LBL_PISA_COMMONKMERS_HPP

#include "../Types.hpp"
namespace pisa{
  struct CommonKmers {
    /*! The number of common kmers between two sequences.
     * The maximum could be floor((l-k)/s)+1, where
     * l is the sequence length, k is the kmer length, and
     * s is the stride. Since l is within 2^16-1 (unsigned short max)
     * we can represent the count as unsigned short as well.
     */
    ushort count;

    /*! The position within the sequence, which is
     * much less than 2^16 - 1 for proteins
     */
    std::pair<ushort, ushort> first;
    std::pair<ushort, ushort> second;

    CommonKmers() : count(1) {
    }

    explicit CommonKmers(ushort count) : count(count){
    }

    friend std::ostream &operator<<(std::ostream &os, const CommonKmers &m) {
      os << "|" << m.count << "(" << m.first.first << "," << m.first.second
         << ")(" <<
         m.second.first << "," << m.second.second << ")| ";
      return os;
    }
  };
}
#endif //LBL_PISA_COMMONKMERS_HPP
