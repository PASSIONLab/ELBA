#ifndef PAIRWISE_ALIGNMENT_H_
#define PAIRWISE_ALIGNMENT_H_

#include "common.h"
#include "DistributedFastaData.hpp"
#include "SharedSeeds.hpp"
#include "Overlap.hpp"

CT<Overlap>::PSpParMat PairwiseAlignment(DistributedFastaData& dfd, CT<SharedSeeds>::PSpParMat& Bmat, int mat, int mis, int gap, int dropoff);

#endif
