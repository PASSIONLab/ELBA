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

// void DistributedFastaData::findnbrs(std::vector<NbrData>& nbrs, int rc_flag)
// {
    // Grid commgrid = index->getcommgrid();
    // int myrank = commgrid->GetRank();
    // int nprocs = commgrid->GetSize();
    // MPI_Comm comm = commgrid->GetWorld();
// }
