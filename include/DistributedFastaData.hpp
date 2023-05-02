#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaIndex.hpp"
#include "Logger.hpp"

class DistributedFastaData
{
public:

    struct FastaDataRequest
    {
        int owner, requester; /* owner is orignal data owner, requester is remote requester */
        size_t offset, count; /* offset is global sequence offset, count is the number of sequences requested */
        unsigned short rc; /* rc=0 means sequences are common cross the row processor dimension, rc=1 ... column processor dimension */

        FastaDataRequest() = default;
        FastaDataRequest(int owner, int requester, size_t offset, size_t count, unsigned short rc)
            : owner(owner), requester(requester), offset(offset), count(count), rc(rc) {}

        friend std::ostream& operator<<(std::ostream& stream, const FastaDataRequest& req)
        {
            stream << "{owner=" << req.owner << ";requester=" << req.requester << ";offset=" << req.offset << ";count=" << req.count << ";rc=" << req.rc << "}";
            return stream;
        }
    };

    DistributedFastaData(std::shared_ptr<FastaIndex> index);

    size_t getrowstartid() const { return rowstartid; }
    size_t getcolstartid() const { return colstartid; }
    size_t getnumrowreads() const { return numrowreads; }
    size_t getnumcolreads() const { return numcolreads; }

    void collect_row_sequences(std::shared_ptr<DnaBuffer> mydna);
    void collect_col_sequences(std::shared_ptr<DnaBuffer> mydna);
    void wait();

    std::shared_ptr<DnaBuffer> getrowbuf() { return rowbuf; }
    std::shared_ptr<DnaBuffer> getcolbuf() { return colbuf; }

private:
    std::shared_ptr<FastaIndex> index;
    bool isdiag;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;

    std::shared_ptr<DnaBuffer> rowbuf, colbuf;

    size_t rowreqbufsize, rowreqnumreads;
    size_t colreqbufsize, colreqnumreads;
    std::unique_ptr<size_t[]> rowreqreadlens;
    std::unique_ptr<size_t[]> colreqreadlens;
    std::unique_ptr<uint8_t[]> rowreqbuf;
    std::unique_ptr<uint8_t[]> colreqbuf;
    std::vector<MPI_Request> rowsendreqs, rowrecvreqs;
    std::vector<MPI_Request> colsendreqs, colrecvreqs;

    void getremoterequests(std::vector<FastaDataRequest>& allrequests, std::vector<FastaDataRequest>& myrequests) const;
    void getgridrequests(std::vector<FastaDataRequest>& myrequests, size_t globalstartid, size_t count, unsigned short rc) const;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
