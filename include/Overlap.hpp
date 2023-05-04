#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <limits>
#include <cassert>
#include "SharedSeeds.hpp"
#include "XDropAligner.hpp"

class Overlap
{
public:
    Overlap() = default;
    Overlap(SeedPair len, SeedPair seed) : beg{}, end{}, len(len), seed(seed), score(0), suffix(0), suffixT(0), direction(-1), directionT(-1) { SetPathInf(); }

    void extend_overlap(const DnaSeq& seqQ, const DnaSeq& seqT, int mat, int mis, int gap, int dropoff);


private:
    SeedPair beg, end, len;
    SeedPair seed;
    int score;
    int suffix, suffixT;
    int8_t direction, directionT;
    int suffix_paths[4];

    void SetPathInf() { std::fill_n(suffix_paths, 4, std::numeric_limits<int>::max()); }
};

void Overlap::extend_overlap(const DnaSeq& seqQ, const DnaSeq& seqT, int mat, int mis, int gap, int dropoff)
{
    assert(seqQ.size() == std::get<0>(len) && seqT.size() == std::get<1>(len));

    XSeed result;

    xdrop_aligner(seqQ, seqT, std::get<0>(seed), std::get<1>(seed), mat, mis, gap, dropoff, result);
}

#endif
