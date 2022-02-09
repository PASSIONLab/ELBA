
#ifndef TRANSITIVE_REDUCTION_H_
#define TRANSITIVE_REDUCTION_H_

#include "../include/ReadOverlap.hpp"

#include <sys/time.h>
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>
#include <map>
#include <fstream>
#include <string>
#include <sstream>


void TransitiveReduction(SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>>& R)
{
    SpParMat<int64_t, ReadOverlap, SpDCCols<int64_t, ReadOverlap>> RT = R;
    RT.Transpose();
    RT.Apply(TransposeSRing());

    // RT.ParallelWriteMM("read_overlaps_symmetric.mm", true, ReadOverlapHandler());

    if (!(RT == R)) R += RT;

    R.Prune(InvalidSRing());

    R.ParallelWriteMM("read_overlaps_symmetric.mm", true, ReadOverlapMMHandler());

}


#endif
