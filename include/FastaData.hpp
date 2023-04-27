#ifndef LBL_DAL_FASTADATA_HPP
#define LBL_DAL_FASTADATA_HPP

#include "common.h"
#include "FastaIndex.hpp"
#include <cassert>


/*!
 * Utility to store and retrive the data read from FASTA files.
 */
class FastaData
{
public:
    FastaData(FIndex index);
    ~FastaData() { delete[] buf; }
    void log() const;

    Grid getcommgrid() const { return index->getcommgrid(); }
    FIndex getindex() const { return index; }
    size_t getfirstid() const { return index->getsomefirstid(idxtag); }
    std::string getsequence(size_t localid) const;
    void ParallelWrite(const char *fname) const;

private:
    FIndex index;
    std::vector<size_t> byteoffsets;
    std::vector<size_t> readlens;
    uint8_t *buf;
    int idxtag; /* corresponds to the processor first responsible for these sequences */

    size_t getbufsize() const { return byteoffsets.back(); }

    static size_t bytesneeded(size_t numbases);
    static size_t computebufbound(size_t totbases, size_t numreads);
};

#endif //LBL_DAL_FASTADATA_HPP
