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

    FIndex getindex() const { return index; }

private:
    FIndex index;
    std::vector<size_t> byteoffsets;
    std::vector<size_t> readlens;
    uint8_t *buf;

    size_t getbufsize() const { return byteoffsets.back(); }

    static size_t bytesneeded(size_t numbases);
    static size_t computebufbound(size_t totbases, size_t numreads);
};

#endif //LBL_DAL_FASTADATA_HPP
