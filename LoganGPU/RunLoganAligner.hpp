/* Created by Giulia Guidi on 4/20/2021. */

#ifndef __LOGAN_ALIGN_RESULT_HPP__
#define __LOGAN_ALIGN_RESULT_HPP__

#include <omp.h>
#include <chrono>
#include <cmath>
#include <string>
#include <vector>

#include "logan.hpp"
#include "interface.hpp"

#define BATCH_SIZE 100000

void RunLoganAlign(vector<string>& seqHs, vector<string>& seqVs, 
	vector<SeedInterface>& seeds, vector<LoganResult>& xscores, int& xdrop, ushort& seed_length);

#endif // __LOGAN_ALIGN_RESULT_HPP__