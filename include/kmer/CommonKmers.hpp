// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef DIBELLA_COMMONKMERS_HPP
#define DIBELLA_COMMONKMERS_HPP

#include "../Types.hpp"
#include "../Defines.hpp"

// for benchmarking
#define EXTRA

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
	 *	StringMatrixEntry econdes both direction and overhang length for both strands */
	// std::vector<uint32_t> overhang(2, 0);
	uint32_t overhang;
	uint32_t overhangT;

#ifdef EXTRA
	uint32_t lenv;
	uint32_t lenh;
	uint32_t overlap;
#endif

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

	operator bool() const { return overhang; };

    // Overload + operator to add two CommonKmers objects
	// Used for: B += BT (TransitiveReductionSR.hpp)
	// The +operator when creating the symmetric matrix doesn't really matter, there's gonna be zeros on the other side
	// The +operator might be needed (meaningful) elsewhere later
    CommonKmers operator+(const CommonKmers& b)
	{
		CommonKmers myobj;
		myobj = b;
		return myobj;
    }

	// Used for: if(!(BT == B)) (TransitiveReductionSR.hpp)
	friend bool operator==(const CommonKmers& lhs, const CommonKmers& rhs)
	{
		if(lhs.overhang == rhs.overhang) return true;
		else return false;
	}

	// Used in SR.hpp fo MinPlus
	// friend bool operator<(const CommonKmers& lhs, const CommonKmers& rhs)
	// {
	// 	ushort len1 = lhs.overhang >> 2;
	// 	ushort len2 = rhs.overhang >> 2;

	// 	if(len1 < len2) return true;
	// 	else return false;
	// }

	// Used in SR.hpp fo MinPlus
	// friend CommonKmers operator+(const CommonKmers& lhs, const CommonKmers& rhs)
	// {
	// 	CommonKmers me;

	// 	ushort dir;

	// 	int mybin1[2] = {0, 0};
	// 	int mybin2[2] = {0, 0};

	// 	if((lhs.overhang & 3) != 0)
	// 	{
	// 		int nbit = 2;
	// 		uint n = lhs.overhang & 3;
	// 		for(int i = 0; i < nbit; i++)
	// 		{
	// 			mybin1[i] = n % 2;
	// 			n = n / 2;
	// 		}
	// 	}

	// 	if((rhs.overhang & 3) != 0)
	// 	{
	// 		int nbit = 2;
	// 		uint n = rhs.overhang & 3;
	// 		for(int i = 0; i < nbit; i++)
	// 		{
	// 			mybin2[i] = n % 2;
	// 			n = n / 2;
	// 		}
	// 	}

	// 	ushort start = mybin1[1];
	// 	ushort end   = mybin2[0];

	// 	if(start == 0)
	// 	{
	// 		if(end == 0) dir = 0;
	// 		else dir = 1;
	// 	}
	// 	else
	// 	{
	// 		if(end == 0) dir = 2;
	// 		else dir = 3;
	// 	}

	// 	ushort len1 = lhs.overhang >> 2;
	// 	ushort len2 = rhs.overhang >> 2;

	// 	len1 += len2;

	// 	me.overhang = len1 << 2 | dir;
	// 	return me;
	// }

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

	struct CkOutputHandler
	{
		template <typename c, typename t>
		void save(std::basic_ostream<c,t> &os,
				const dibella::CommonKmers &v,
				uint64_t row,
				uint64_t col)
		{
			// GGGG: we need the overhand value to create input in graph dot for comparison
			int dir = v.overhang & 3;
			int rc  = 0;
			if(dir == 0 || dir == 3) rc = 1;

			// direction, rc, overhang, begV, endV, begH, endH (OverlapLen and others computed in python script during translation)
			os << dir << "\t" << rc << "\t" << v.overhang << "\t" << v.first.first << "\t" << v.first.second << "\t"
				<< v.second.first << "\t" <<
				#ifdef EXTRA
				v.second.second  << "\t"  << v.lenv << "\t" << v.lenh << "\t" << v.overlap;
				#else
				v.second.second;
				#endif
		}
	};

    struct CkOutputMMHandler
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t> &os,
                        const dibella::CommonKmers &v,
                        uint64_t row,
                        uint64_t col)
        {
			int dir = v.overhang  & 3;
			int len = v.overhang >> 2;
			os << dir << "\t" << len;
			// os;
			// std::string val = std::to_string(len) + "\t" + std::to_string(dir);
			// os << val;
        }
    };

    struct CkOutputMMHandlerBool
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t> &os,
                        const bool &v,
                        uint64_t row,
                        uint64_t col)
        {
            os << v;
        }
    };

}

#endif //DIBELLA_COMMONKMERS_HPP
