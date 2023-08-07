#ifndef XDROP_ALIGNER_H_
#define XDROP_ALIGNER_H_

#include "common.h"
#include "DnaSeq.hpp"
#include <vector>
#include <algorithm>

#define DELTACHERNOFF (0.1)

typedef enum
{
    BAD_ALIGNMENT,
    FIRST_CONTAINED,
    SECOND_CONTAINED,
    FIRST_TO_SECOND_OVERLAP,
    SECOND_TO_FIRST_OVERLAP
} OverlapClass;

struct XSeed
{
    int begQ, endQ, begT, endT, score;
    bool rc;

    XSeed() : begQ(0), endQ(0), begT(0), endT(0), score(-1), rc(false) {}
};

int xdrop_aligner(const DnaSeq& seqQ, const DnaSeq& seqT, int begQ, int begT, int mat, int mis, int gap, int dropoff, XSeed& result);
void classify_alignment(const XSeed& ai, int lenQ, int lenT, float target_identity, OverlapClass& kind);

#endif
