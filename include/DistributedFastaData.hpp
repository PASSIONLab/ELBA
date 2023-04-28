#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"

class DistributedFastaData
{
public:
    struct NbrData { uint32_t id : 31, rc_flag : 1; };
    static_assert(sizeof(NbrData) == 4);

    static NbrData nbrdata_ctor(int rank, int rc_flag) { return {static_cast<uint32_t>(rank), static_cast<uint32_t>(!!rc_flag)}; }

    DistributedFastaData(FIndex index);

    void findnbrs(std::vector<NbrData>& nbrs);

private:
    FIndex index;
    int myrowid, mycolid, procdim;
    size_t readsperprocdim;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
