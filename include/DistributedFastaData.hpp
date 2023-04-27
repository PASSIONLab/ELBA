#ifndef LBL_DAL_DISTRIBUTEDFASTADATA_H
#define LBL_DAL_DISTRIBUTEDFASTADATA_H

#include "FastaData.hpp"
#include "Logger.hpp"

class DistributedFastaData
{
public:

    struct NbrData
    {
        size_t startid;
        size_t numseqs;
        int source;
        int rc_flag;

        NbrData(size_t startid, size_t numseqs, int source, int rc_flag) : startid(startid), numseqs(numseqs), source(source), rc_flag(rc_flag) {}

        // static size_t findnbrs(std::vector<NbrData>& neighbors, const std::vector<size_t>& idoffsets, size_t startid, size_t totreads, int rc_flag, Grid commgrid)
        // {
            // int myrank = commgrid->GetRank();
            // size_t numadded = 0;

            // auto iditr = std::upper_bound(idoffsets.cbegin(), idoffsets.cend(), startid);
            // iditr--;

            // assert(*iditr <= startid);

            // while (*iditr < startid + totreads)
            // {
                // int reqrank = iditr - idoffsets.cbegin();
                // size_t reqstart = std::max(*iditr++, startid);
                // size_t reqend = std::min(*iditr, startid + totreads);
                // neighbors.emplace_back(reqstart, reqend-reqstart, reqrank, rc_flag);
                // numadded++;
            // }

            // return numadded;
        // }
    };

    DistributedFastaData(FIndex index);
    ~DistributedFastaData();

private:
    FIndex index;
    int myrowid, mycolid, procdim;
    size_t readsperprocdim;
    size_t rowstartid, colstartid;
    size_t numrowreads, numcolreads;
};

#endif //LBL_DAL_DISTRIBUTEDFASTADATA_H
