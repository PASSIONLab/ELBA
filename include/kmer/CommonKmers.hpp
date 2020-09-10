// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef DIBELLA_COMMONKMERS_HPP
#define DIBELLA_COMMONKMERS_HPP

#include "../Types.hpp"
#include "../Defines.hpp"

namespace dibella {
  struct CommonKmers {
    /*! The number of common kmers between two sequences.
     * The maximum could be floor((l-k)/s)+1, where
     * l is the sequence length, k is the kmer length, and
     * s is the stride. Since l is within 2^16-1 (unsigned short max)
     * we can represent the count as unsigned short as well.
     */
    ushort count;

	bool passed;
	uint32_t score; /* Used for storing alignment score */
	
	/*! GGGG: this is either the suffix or prefix entry need for the transitive reduction 
	 *	StringMatrixEntry econdes both direction and overhang length*/
	uint32_t overhang; 

    /*! The position within the sequence, which is
     * much less than 2^16 - 1 for proteins
     */

#ifdef TWOSEED
	// GGGG: just use two seeds per read
    std::pair<PosInRead, PosInRead> first;
    std::pair<PosInRead, PosInRead> second;
#else
	// GGGG: need this to compute distance
	std::vector<std::pair<PosInRead, PosInRead>> pos;
#endif

    CommonKmers() : count(1), passed(false), overhang(0) {
    }
    explicit
	CommonKmers(ushort count) : 
		count(count), passed(false), overhang(0) {
    }

	CommonKmers (bool passed, uint32_t score) :
		passed(passed),
		score(score) {
	}

    friend std::ostream &operator<<(std::ostream &os, const CommonKmers &m)
	{
	#ifdef TWOSEED
		os << "|" << m.count << "(" << m.first.first << "," << m.first.second
			<< ")(" <<
			m.second.first << "," << m.second.second << ")| ";
	#else
		os << "|" << m.count << "(";
		for(int i = 0; i < m.pos.size(); i++)
		{
			os << m.pos[i].first << "," << m.pos[i].second << ")| ";  
		}
	#endif
		return os;
    }

	};

	/* GGGG: matrix symmetrication removed */

	struct CkOutputHandler
	{
		template <typename c, typename t>
		void save(std::basic_ostream<c,t> &os,
				const dibella::CommonKmers &v,
				uint64_t row,
				uint64_t col)
		{
			os << v.score;
		}
	};
}

#endif //DIBELLA_COMMONKMERS_HPP
