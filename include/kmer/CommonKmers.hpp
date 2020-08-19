// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef DIBELLA_COMMONKMERS_HPP
#define DIBELLA_COMMONKMERS_HPP

#include "../Types.hpp"

namespace dibella {
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

#ifdef TWOSEED
	// GGGG: just use two seeds per read
    std::pair<ushort, ushort> first;
    std::pair<ushort, ushort> second;
#else
	// GGGG: need this to compute distance
	std::vector<std::pair<PosInRead, PosInRead>> pos;
#endif

    CommonKmers() : count(1) {
    }

    explicit CommonKmers(ushort count) : count(count){
    }

#ifndef TWOSEED
	explicit CommonKmers(CommonKmers kmer) : count(kmer.count), pos(kmer.pos) {
    }
#endif

    friend std::ostream &operator<<(std::ostream &os, const CommonKmers &m)
	{
#ifdef TWOSEED
      os << "|" << m.count << "(" << m.first.first << "," << m.first.second
         << ")(" <<
         m.second.first << "," << m.second.second << ")| ";
#else
      os << "|" << m.count << "(";
	  for(int i = 0; i < pos.size(); i++)
	  {
	  	os << m.pos[i].first << "," << m.pos[i].second << ")| ";  
	  }
#endif
      return os;
    }

	// used in the symmetricization of the output matrix
	CommonKmers
	operator+(const CommonKmers &arg)
	{
		CommonKmers tmp(0);
	#ifdef TWOSEED
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
		else // both should have count = 1
		{
			tmp.count  = 2;
			tmp.first  = this->first;
			tmp.second = arg.first;
		}
	#else
		if (this->count >= 2)
		{
			tmp.count = this->count + arg.count;
			tmp.pos   = this->pos;
		}
		else if (arg.count >= 2)
		{
			tmp.count = this->count + arg.count;
			tmp.pos   = arg.pos;
		}
		else // both should have count = 1
		{
			tmp.count = 2;
			tmp.pos   = this->pos;
		}
	#endif
		return tmp;
	}
  };
}
#endif //DIBELLA_COMMONKMERS_HPP
