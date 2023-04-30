#ifndef LBL_DAL_FASTADATA_HPP
#define LBL_DAL_FASTADATA_HPP

#include "common.h"
#include "FastaIndex.hpp"
#include "DnaSeq.hpp"
#include <cassert>


/*!
 * Utility to store and retrive the data read from FASTA files.
 */
class FastaData
{
public:
    FastaData(FIndex index);
    ~FastaData() = default;
    void log(FIndex index) const;

private:
    std::unique_ptr<DnaBuffer> buffer;
    std::vector<DnaSeq> sequences;
};

#endif //LBL_DAL_FASTADATA_HPP
