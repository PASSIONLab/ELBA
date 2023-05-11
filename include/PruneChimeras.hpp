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

    PileupVector(int read_length);
    PileupVector(const std::vector<int>& v, int offset, int size);

    int Length() const;

    void AddInterval(int begin, int end);

    ~PileupVector() = default;
};

std::vector<int> PackPileupVector(const std::vector<PileupVector>& pvs, std::vector<int>& lens, int& size);
std::vector<PileupVector> UnpackPileupVector(const std::vector<int>& packed, const std::vector<int>& lens);
/* void add_gaps(int begin, int end, int length, std::vector<std::tuple<int, int>>& middle, std::vector<std::tuple<int, int>>& extremity); */

std::vector<PileupVector> GetReadPileup(DistributedFastaData& dfd, CT<Overlap>::PSpParMat& Rmat);
/* FullyDistVec<int64_t, int64_t> GetChimeras(const std::shared_ptr<CommGrid>& comgrid, std::shared_ptr<DistributedFastaData> dfd, const std::vector<PileupVector>& pileups, int coverage_min); */

#endif
