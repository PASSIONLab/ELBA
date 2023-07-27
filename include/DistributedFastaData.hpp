#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaIndex.hpp"
#include "DnaBuffer.hpp"
#include <iostream>
#include <memory>
#include <cstddef>

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
            stream << "{" << req.owner+1 << "," << req.requester+1 << "," << req.offset << "," << req.count << "," << req.rc << "}";
            return stream;
        }
    };

    DistributedFastaData(FastaIndex& index);

    FastaIndex& getindex() { return index; }

    size_t getrowstartid() const { return rowinfo.startid; }
    size_t getcolstartid() const { return colinfo.startid; }
    size_t getnumrowreads() const { return rowinfo.numreads; }
    size_t getnumcolreads() const { return colinfo.numreads; }

    void collect_sequences(const DnaBuffer& mydna);
    void wait();

    std::shared_ptr<DnaBuffer> getrowbuf() { return rowbuf; }
    std::shared_ptr<DnaBuffer> getcolbuf() { return colbuf; }

    void write_grid_sequences(char const *fname_prefix) const;

    struct DimExchangeInfo
    {
        size_t startid, numreads;
        size_t reqbufsize, reqnumreads;
        std::unique_ptr<size_t[]> reqreadlens;
        std::unique_ptr<uint8_t[]> reqbuf;
        std::vector<MPI_Request> sendreqs, recvreqs;
    };

private:
    FastaIndex& index;
    bool isdiag;

    std::shared_ptr<DnaBuffer> rowbuf, colbuf;
    DimExchangeInfo rowinfo, colinfo;

    void collect_dim_sequences(const DnaBuffer& mydna, DimExchangeInfo& diminfo);
    void getremoterequests(std::vector<FastaDataRequest>& allrequests, std::vector<FastaDataRequest>& myrequests) const;
    void getgridrequests(std::vector<FastaDataRequest>& myrequests, size_t globalstartid, size_t count, unsigned short rc) const;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
