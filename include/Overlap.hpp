#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <limits>
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
};

#endif
