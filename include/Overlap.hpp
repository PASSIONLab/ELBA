#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <limits>
#include <iostream>
#include <cassert>
#include "SharedSeeds.hpp"

struct Overlap
{
    Overlap(SeedPair len, SeedPair seed) :
        beg{}, end{}, len(len),
        seed(seed),
        score(0),
        suffix(0), suffixT(0),
        direction(-1), directionT(-1),
        rc(false), passed(false), containedQ(false), containedT(false) { SetPathInf(); }

    Overlap(const Overlap& rhs) :
        beg(rhs.beg), end(rhs.end), len(rhs.len),
        seed(rhs.seed),
        score(rhs.score),
        suffix(rhs.suffix), suffixT(rhs.suffixT),
        direction(rhs.direction), directionT(rhs.directionT),
        rc(rhs.rc), passed(rhs.passed), containedQ(rhs.containedQ), containedT(rhs.containedT) { std::copy(rhs.suffix_paths, rhs.suffix_paths + 4, suffix_paths); }

    Overlap() : Overlap({}, {}) {}

    void extend_overlap(const DnaSeq& seqQ, const DnaSeq& seqT, int mat, int mis, int gap, int dropoff);
    void classify();

    SeedPair beg, end, len;
    SeedPair seed;
    int score;
    int suffix, suffixT;
    int8_t direction, directionT;
    int suffix_paths[4];
    bool rc, passed, containedQ, containedT;

    void SetPathInf() { std::fill_n(suffix_paths, 4, std::numeric_limits<int>::max()); }

    friend std::ostream& operator<<(std::ostream& os, const Overlap& o)
    {
        char rcflag = o.rc? '-' : '+';
        os << std::get<0>(o.len) << "\t" << std::get<0>(o.beg) << "\t" << std::get<0>(o.end) << "\t" << rcflag << "\t" << std::get<1>(o.len) << "\t" << std::get<1>(o.beg) << "\t" << std::get<1>(o.end) << "\t" << o.score;
        return os;
    }

    struct Transpose : std::unary_function<Overlap, Overlap>
    {
        Overlap operator()(const Overlap& o) const
        {
            Overlap oT(o);
            PosInRead tmp;

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
};

#endif
