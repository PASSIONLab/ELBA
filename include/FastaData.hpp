// Created by Saliya Ekanayake on 1/7/19.

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
    FastaData(Grid commgrid) : commgrid(commgrid), byteoffsets(1, 0) {}
    FastaData(const FastaIndex& index);
    FastaData(const FastaData& rhs);

    size_t numreads() const { assert(readlens.size() == byteoffsets.size()-1); return readlens.size(); }
    size_t getfirstid() const { return firstid; }
    Grid getcommgrid() const { return commgrid; }
    size_t getbufsize() const { return byteoffsets.back(); }
    void log() const;

    std::string getsequence(size_t localid) const;
    uint8_t const* getbufrange(size_t localstartid, size_t rangelen, size_t& numbytes) const;

    ~FastaData() { delete[] buf; }

    static size_t bytesneeded(size_t numbases);
    static size_t computebufbound(size_t totbases, size_t numreads);

private:
    Grid commgrid;
    uint8_t *buf;                    /* buffer storing local read sequences (compressed) */
    std::vector<size_t> byteoffsets; /* offsets within compressed buf of first byte in a read sequence */
    std::vector<size_t> readlens;    /* uncompressed read lengths */
    size_t firstid;                  /* these represent global read ids */
};


#endif //LBL_DAL_FASTADATA_HPP
