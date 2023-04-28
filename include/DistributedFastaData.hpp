#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"
#include "DnaSeq.hpp"

class DistributedFastaData
{
public:
    struct NbrData
    {
        size_t startid, numseqs;
        int sendrank, destrank;
        unsigned short rcflag;

        NbrData() = default;
    };

    DistributedFastaData(FIndex index);

    void findnbrs(std::vector<NbrData>& mynbrs, size_t startid, size_t count, unsigned short rc_flag) const;

    void allgather_neighbors();
    void exchange_reads();

private:
    FIndex index;
    int myrowid, mycolid, procdim;
    size_t readsperprocdim;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;

    uint8_t *rowbuf, *colbuf;
    std::vector<DnaSeq> rowreads, colreads;
    std::vector<NbrData> allnbrs;
};

typedef typename DistributedFastaData::NbrData NbrData;

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
