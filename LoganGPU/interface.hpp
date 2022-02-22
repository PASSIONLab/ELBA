#ifndef __LOGAN_INTERFACE_CUH__
#define __LOGAN_INTERFACE_CUH__

#include <omp.h>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

/*
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
*/
#endif
