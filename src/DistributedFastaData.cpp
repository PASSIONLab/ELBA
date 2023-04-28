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
    readsperprocdim = numreads / procdim;
    rowstartid = myrowid * readsperprocdim;
    colstartid = mycolid * readsperprocdim;

    numrowreads = (myrowid == procdim-1)? (numreads - myrowid*readsperprocdim) : readsperprocdim;
    numcolreads = (mycolid == procdim-1)? (numreads - mycolid*readsperprocdim) : readsperprocdim;
}

void DistributedFastaData::findnbrs(std::vector<NbrData>& nbrs)
{
    Grid commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    std::vector<int> owners;

    owners = index->collectowners(rowstartid, numrowreads);

    for (auto itr = owners.begin(); itr != owners.end(); ++itr)
    {
        nbrs.push_back(nbrdata_ctor(*itr, 0));
    }

    if (myrowid != mycolid)
    {
        owners = index->collectowners(colstartid, numcolreads);

        for (auto itr = owners.begin(); itr != owners.end(); ++itr)
        {
            nbrs.push_back(nbrdata_ctor(*itr, 1));
        }
    }

    std::sort(nbrs.begin(), nbrs.end(), [](const auto& a, const auto& b) { return a.id < b.id; });

    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << "row range " << Logger::readrangestr(rowstartid, numrowreads) << ". "
                                                   << "col range " << Logger::readrangestr(colstartid, numcolreads) << "\n";

    for (auto itr = nbrs.begin(); itr != nbrs.end(); ++itr)
    {
        logger() << " requests sequences from " << logger.rankstr(itr->id) << ": " << Logger::readrangestr(index->getsomefirstid(itr->id), index->getreadcount(itr->id)) << ". (rc_flag=" << itr->rc_flag << ")\n";
    }
    logger.Flush("Exchange parameters:");
}
