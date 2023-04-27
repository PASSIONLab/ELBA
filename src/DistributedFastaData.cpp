#include "DistributedFastaData.hpp"
#include "Logger.hpp"

DistributedFastaData::DistributedFastaData(FIndex index) : index(index)
{
    Grid commgrid = index->getcommgrid();

    assert(commgrid->GetGridRows() == commgrid->GetGridCols());

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    myrowid = commgrid->GetRankInProcCol();
    mycolid = commgrid->GetRankInProcRow();
    procdim = commgrid->GetGridRows();

    size_t numreads = index->gettotrecords();
    rowstartid = myrowid * readsperprocdim;
    colstartid = mycolid * readsperprocdim;
    numrowreads = (myrowid == procdim-1)? (numreads - myrowid*readsperprocdim) : readsperprocdim;
    numcolreads = (mycolid == procdim-1)? (numreads - mycolid*readsperprocdim) : readsperprocdim;
}

    // struct NbrData
    // {
        // size_t startid;
        // size_t numseqs;
        // int source;
        // int rc_flag;

        // NbrData(size_t startid, size_t numseqs, int source, int rc_flag) : startid(startid), numseqs(numseqs), source(source), rc_flag(rc_flag) {}

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
    // };
