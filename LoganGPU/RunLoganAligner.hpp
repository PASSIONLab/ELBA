/* Created by Giulia Guidi on 4/20/2021. */

#ifndef __LOGAN_ALIGN_RESULT_HPP__
#define __LOGAN_ALIGN_RESULT_HPP__

#include <omp.h>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#define BATCH_SIZE 100000

using namespace std;
struct LoganResult {

    int  score;
    bool    rc;

    int begSeedH;
    int begSeedV;
    int endSeedH;
    int endSeedV;
};

struct SeedInterface {

        int beginPositionH;
        int beginPositionV;

        int endPositionH;
        int endPositionV;

        int seedLength;
        int score;

        SeedInterface(): beginPositionH(0), beginPositionV(0), endPositionH(0), endPositionV(0), score(0)
        { }

        SeedInterface(int beginPositionH, int beginPositionV, int seedLength):
                beginPositionH(beginPositionH), beginPositionV(beginPositionV), endPositionH(beginPositionH + seedLength),
                endPositionV(beginPositionV + seedLength), score(0)
        { }

        SeedInterface(int beginPositionH, int beginPositionV, int endPositionH, int endPositionV):
                beginPositionH(beginPositionH),
                beginPositionV(beginPositionV),
                endPositionH(endPositionH),
                endPositionV(endPositionV),
                score(0)
        { }

        SeedInterface(SeedInterface const& other):
                beginPositionH(other.beginPositionH),
                beginPositionV(other.beginPositionV),
                endPositionH(other.endPositionH),
                endPositionV(other.endPositionV),
                score(0)
        { }
};

void 
RunLoganAlign(vector<string>& seqHs, vector<string>& seqVs, 
	vector<SeedInterface>& SeedInterfaceSet, vector<LoganResult>& xscores, int& xdrop, ushort& seed_length);

#endif // __LOGAN_ALIGN_RESULT_HPP__