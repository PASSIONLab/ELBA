#include "SharedSeeds.hpp"
#include "common.h"

std::unique_ptr<CT<SharedSeeds>::PSpParMat>
create_seed_matrix(CT<PosInRead>::PSpParMat& A, CT<PosInRead>::PSpParMat& AT)
{
    auto B = std::make_unique<CT<SharedSeeds>::PSpParMat>(Mult_AnXBn_DoubleBuff<SharedSeeds::Semiring, SharedSeeds, CT<SharedSeeds>::PSpDCCols>(A, AT));
    B->Prune([](const SharedSeeds& nz) { return nz.getnumshared() <= 1; });
    return std::move(B);
}
