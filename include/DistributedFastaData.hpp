#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"

class DistributedFastaData
{
public:

    struct FastaDataRequest
    {
        int owner, requester; /* owner is orignal data owner, requester is remote requester */
        size_t offset, count; /* offset is local to the data owned by the requester, count is the number of sequences requested */
        unsigned short rc; /* rc=0 means sequences are common cross the row processor dimension, rc=1 ... column processor dimension */

        FastaDataRequest(int owner, int requester, size_t offset, size_t count, unsigned short rc)
            : owner(owner), requester(requester), offset(offset), count(count), rc(rc) {}
    };

    DistributedFastaData(std::shared_ptr<FastaIndex> index);

    size_t getrowstartid() const { return rowstartid; }
    size_t getcolstartid() const { return colstartid; }
    size_t getnumrowreads() const { return numrowreads; }
    size_t getnumcolreads() const { return numcolreads; }

private:
    std::shared_ptr<FastaIndex> index;
    bool isdiag;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;

    std::vector<FastaDataRequest> getremoterequests() const;
    void getgridrequests(std::vector<FastaDataRequest>& myrequests, size_t globalstartid, size_t count, unsigned short rc) const;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
