#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"

class DistributedFastaData
{
public:

    struct FastaDataRequest
    {
        int sender, receiver; /* sender is orignal data owner, receiver is remote requester */
        size_t offset, count; /* offset is local to the data owned by the receiver, count is the number of sequences requested */
        unsigned short rc; /* rc=0 means sequences are common cross the row processor dimension, rc=1 ... column processor dimension */

        FastaDataRequest(int sender, int receiver, size_t offset, size_t count, unsigned short rc)
            : sender(sender), receiver(receiver), offset(offset), count(count), rc(rc) {}
    };

    DistributedFastaData(FIndex index);

    size_t getrowstartid() const { return rowstartid; }
    size_t getcolstartid() const { return colstartid; }
    size_t getnumrowreads() const { return numrowreads; }
    size_t getnumcolreads() const { return numcolreads; }

private:
    FIndex index;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
