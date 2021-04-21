/* Created by Giulia Guidi on 4/20/2021. */

#ifndef __LOGAN_ALIGN_RESULT_HPP__
#define __LOGAN_ALIGN_RESULT_HPP__

#include <omp.h>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#define BATCH_SIZE 100000

struct LoganResult {
    int  score;
    bool    rc;

    int begSeedH;
    int begSeedV;
    int endSeedH;
    int endSeedV;
};

#endif // __LOGAN_ALIGN_RESULT_HPP__