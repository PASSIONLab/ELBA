#include "Overlap.hpp"
#include "XDropAligner.hpp"

Overlap::Overlap(std::tuple<PosInRead, PosInRead> len, std::tuple<PosInRead, PosInRead> seed) :
    beg{}, end{}, len(len),
    seed(seed),
    score(0),
    suffix(0), suffixT(0),
    direction(-1), directionT(-1),
    rc(false), passed(false), containedQ(false), containedT(false) { SetPathInf(); }

Overlap::Overlap(const Overlap& rhs) :
    beg(rhs.beg), end(rhs.end), len(rhs.len),
    seed(rhs.seed),
    score(rhs.score),
    suffix(rhs.suffix), suffixT(rhs.suffixT),
    direction(rhs.direction), directionT(rhs.directionT),
    rc(rhs.rc), passed(rhs.passed), containedQ(rhs.containedQ), containedT(rhs.containedT) { std::copy(rhs.suffix_paths, rhs.suffix_paths + 4, suffix_paths); }

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

    int lenQ  = std::get<0>(len);
    int lenT  = std::get<1>(len);
    int begQr = result.begQ;
    int endQr = result.endQ;
    int begTr = rc? std::get<1>(len) - result.endT : result.begT;
    int endTr = rc? std::get<1>(len) - result.begT : result.endT;

    if (kind != BAD_ALIGNMENT)
    {
        passed = true;

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
            suffix     = ((lenT - endTr) - (lenQ - endQr));
            suffixT    = begQr - begTr;
        }
        else
        {
            direction  = rc? 3 : 2;
            directionT = rc? 3 : 1;
            suffix     = begTr - begQr;
            suffixT    = (lenQ - endQr) - (lenT - endTr);
        }
    }
}
