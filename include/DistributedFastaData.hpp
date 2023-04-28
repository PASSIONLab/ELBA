#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"
#include "DnaSeq.hpp"

class DistributedFastaData
{
public:
    struct NbrData { uint32_t id : 31, rc_flag : 1; };
    static_assert(sizeof(NbrData) == 4);

    static NbrData nbrdata_ctor(int rank, int rc_flag) { return {static_cast<uint32_t>(rank), static_cast<uint32_t>(!!rc_flag)}; }

    DistributedFastaData(FIndex index);

    void allgather_neighbors() {}
    void exchange_reads() {}

private:
    FIndex index;
    int myrowid, mycolid, procdim;
    size_t readsperprocdim;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;

    uint8_t *rowbuf, *colbuf;
    std::vector<DnaSeq> rowreads, colreads;
    std::vector<NbrData> allnbrs;

    void findnbrs(std::vector<NbrData>& nbrs);
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
