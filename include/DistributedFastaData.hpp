#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"

class DistributedFastaData
{
public:
    struct NbrData { uint16_t id : 15, rc_flag : 1; };
    static_assert(sizeof(NbrData) == 2);

    DistributedFastaData(FIndex index);
    ~DistributedFastaData();

private:
    FIndex index;
    int myrowid, mycolid, procdim;
    size_t readsperprocdim;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
