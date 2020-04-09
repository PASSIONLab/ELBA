// Created by Saliya Ekanayake on 10/15/19.

#ifndef DISTAL_COMMONKMERS_HPP
#define DISTAL_COMMONKMERS_HPP

#include "../Types.hpp"
namespace distal{
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

	// used in the symmetricization of the output matrix
	CommonKmers
	operator+(const CommonKmers &arg)
	{
		CommonKmers tmp(0);
		if (this->count >= 2)
		{
			tmp.count  = this->count + arg.count;
			tmp.first  = this->first;
			tmp.second = this->second;
		}
		else if (arg.count >= 2)
		{
			tmp.count  = this->count + arg.count;
			tmp.first  = arg.first;
			tmp.second = arg.second;
		}
		else					// both should have count = 1
		{
			tmp.count  = 2;
			tmp.first  = this->first;
			tmp.second = arg.first;
		}
		return tmp;
	}
  };
}
#endif //DISTAL_COMMONKMERS_HPP
