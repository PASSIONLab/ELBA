// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef ELBA_COMMONKMERS_HPP
#define ELBA_COMMONKMERS_HPP

#include "../Types.hpp"
#include "../Defines.hpp"

// for benchmarking
#define EXTRA

namespace elba {
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

    int dir, dirT, sfx, sfxT;

    bool rc;
	uint32_t lenv;
	uint32_t lenh;

#ifdef EXTRA
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

    CommonKmers() : count(1), passed(false), dir(-1) {}

    explicit
	CommonKmers(ushort count) : count(count), passed(false), dir(-1) {}

	CommonKmers (bool passed, uint32_t score) : passed(passed), score(score) {}

	operator bool() const { return (dir != -1); }

    bool is_invalid() const { return (dir == -1); }

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
		if(lhs.dir == rhs.dir && lhs.sfx == rhs.sfx) return true;
		else return false;
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

	struct CkOutputHandler
	{
		template <typename c, typename t>
		void save(std::basic_ostream<c,t> &os,
				const elba::CommonKmers &v,
				uint64_t row,
				uint64_t col)
		{

            os << v.dir << "\t" << static_cast<int>(v.rc) << "\t" << v.first.first << "\t" << v.first.second << "\t" << v.second.first << "\t" << v.second.second << "\t"
               << v.lenv << "\t" << v.lenh << "\t" << v.overlap;

		}
	};

    struct CkOutputMMHandler
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t> &os,
                        const elba::CommonKmers &v,
                        uint64_t row,
                        uint64_t col)
        {	
			// In KmerIntersectSR.hpp we have (where res == cks):
			// 	res.first.first 	// Kmer 1 on argA
			// 	res.first.second 	// Kmer 1 on argB
			// 	res.second.first 	// Kmer 2 on argA
			// 	res.second.second 	// Kmer 2 on argB

			os << v.dir << "\t" << v.first.first << "\t" << v.first.second << "\t" << v.second.first << "\t" << v.second.second;
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

struct CommonKmersGraphHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const elba::CommonKmers& v, int64_t row, int64_t col)
    {
        os << v.score << "\t" << v.lenv << "\t" << v.lenh << "\t" << v.rc;
    }
};

#endif //ELBA_COMMONKMERS_HPP