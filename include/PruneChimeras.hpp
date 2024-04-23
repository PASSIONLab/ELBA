#ifndef PRUNE_CHIMERAS_H_
#define PRUNE_CHIMERAS_H_

#include <vector>
#include <tuple>
#include <memory>
#include "common.h"
#include "Overlap.hpp"
#include "DistributedFastaData.hpp"

class PileupVector
{
public:
    std::vector<int> pileup;
    int num_intervals = 0;

    PileupVector(int read_length);
    PileupVector(const std::vector<int>& v, int offset, int size);

    int Length() const;
    int NumIntervals() const;
    int TallestPoint();

    void AddInterval(int begin, int end);

    std::tuple<int, int> GetTrimmedInterval(int threshold);

    ~PileupVector() = default;
};

std::vector<int> PackPileupVector(const std::vector<PileupVector>& pvs, std::vector<int>& lens, int& size);
std::vector<PileupVector> UnpackPileupVector(const std::vector<int>& packed, const std::vector<int>& lens);
std::vector<PileupVector> GetReadPileup(DistributedFastaData& dfd, CT<Overlap>::PSpParMat& Rmat);

#endif
