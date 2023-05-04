#ifndef OVERLAP_H_
#define OVERLAP_H_

#include <algorithm>
#include <limits>
#include <cassert>
#include "SharedSeeds.hpp"
#include "XDropAligner.hpp"

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

void Overlap::extend_overlap(const DnaSeq& seqQ, const DnaSeq& seqT, int mat, int mis, int gap, int dropoff)
{
    assert(seqQ.size() == std::get<0>(len) && seqT.size() == std::get<1>(len));

    XSeed result;

    xdrop_aligner(seqQ, seqT, std::get<0>(seed), std::get<1>(seed), mat, mis, gap, dropoff, result);

    OverlapClass kind;
    classify_alignment(result, seqQ.size(), seqT.size(), kind);

    rc = result.rc;
    score = result.score;

    std::get<0>(beg) = result.begQ;
    std::get<1>(beg) = result.begT;
    std::get<0>(end) = result.endQ;
    std::get<1>(end) = result.endT;

    int begQr = result.begQ;
    int endQr = result.endQ;
    int begTr = rc? std::get<1>(len) - result.endT : result.begT;
    int endTr = rc? std::get<1>(len) - result.begT : result.endT;

    if (kind != BAD_ALIGNMENT)
    {
        if (kind == FIRST_CONTAINED)
        {
            containedQ = true;
        }
        else if (kind == SECOND_CONTAINED)
        {
            containedT = true;
        }
        else if (kind == FIRST_TO_SECOND_OVERLAP)
        {
            direction  = rc? 0 : 1;
            directionT = rc? 0 : 2;
            suffix     = ((std::get<1>(len) - endTr) - (std::get<0>(len) - endQr));
            suffixT    = begQr - begTr;
            passed     = true;
        }
        else
        {
            direction  = rc? 3 : 2;
            directionT = rc? 3 : 1;
            suffix     = begTr - begQr;
            suffixT    = (std::get<0>(len) - endQr) - (std::get<1>(len) - endTr);
            passed     = true;
        }
    }
}

#endif
