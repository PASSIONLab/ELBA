#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <limits>
#include <iostream>
#include <cassert>
#include "SharedSeeds.hpp"

struct Overlap
{
    Overlap() : Overlap({}, {}) {}
    Overlap(std::tuple<PosInRead, PosInRead> len, std::tuple<PosInRead, PosInRead> seed);
    Overlap(const Overlap& rhs);

    operator int() const { return 1; } /* for creating integer matrix with same nonzero pattern */

    void extend_overlap(const DnaSeq& seqQ, const DnaSeq& seqT, int mat, int mis, int gap, int dropoff);
    void classify();

    std::tuple<PosInRead, PosInRead> beg, end, len;
    std::tuple<PosInRead, PosInRead> seed;
    int score;
    int suffix, suffixT;
    int8_t direction, directionT;
    int suffix_paths[4];
    bool rc, passed, containedQ, containedT;

    void SetPathInf() { std::fill_n(suffix_paths, 4, std::numeric_limits<int>::max()); }

    bool arrows(int& t, int& h) const
    {
        if (direction == -1) return false;

        t = (direction >> 1) & 1;
        h = direction & 1;

        return true;
    }

    struct Transpose : std::unary_function<Overlap, Overlap>
    {
        Overlap operator()(const Overlap& o) const
        {
            Overlap oT(o);

            std::get<0>(oT.beg) = std::get<1>(o.beg);
            std::get<1>(oT.beg) = std::get<0>(o.beg);

            std::get<0>(oT.end) = std::get<1>(o.end);
            std::get<1>(oT.end) = std::get<0>(o.end);

            std::get<0>(oT.len) = std::get<1>(o.len);
            std::get<1>(oT.len) = std::get<0>(o.len);

            oT.suffix = o.suffixT;
            oT.suffixT = o.suffix;

            oT.direction = o.directionT;
            oT.directionT = o.direction;

            oT.containedQ = o.containedT;
            oT.containedT = o.containedQ;

            return oT;
        }
    };

    struct IOHandler
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t>& os, const Overlap& o, uint64_t row, uint64_t col) { os << o; }
    };

    Overlap& operator+=(const Overlap& lhs) { return *this; }
    friend Overlap operator+(const Overlap& lhs, const Overlap& rhs) { return lhs; }
    friend bool operator<(const Overlap& lhs, const Overlap& rhs) { return true; }

    friend std::ostream& operator<<(std::ostream& os, const Overlap& o)
    {
        char rcflag = o.rc? '-' : '+';
        os << std::get<0>(o.len) << "\t" << std::get<0>(o.beg) << "\t" << std::get<0>(o.end) << "\t" << rcflag << "\t" << std::get<1>(o.len) << "\t" << std::get<1>(o.beg) << "\t" << std::get<1>(o.end) << "\t" << o.score << "\t" << static_cast<int>(o.containedQ) << "\t" << static_cast<int>(o.containedT);
        return os;
    }
};

#endif
