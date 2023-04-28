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

void DistributedFastaData::findnbrs(std::vector<NbrData>& mynbrs, size_t startid, size_t count, unsigned short rc_flag) const
{
    Grid commgrid = index->getcommgrid();
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    size_t totreads = index->gettotrecords();
    const auto& readdispls = index->getreaddispls();

    assert(count >= 1 && 0 <= startid && startid + count <= totreads);

    std::vector<int> nbrs;
    auto iditr = std::upper_bound(readdispls.cbegin(), readdispls.cend(), static_cast<MPI_Displ_type>(startid));
    iditr--;

    assert(*iditr <= startid);

    while (*iditr < startid + count)
    {
        int sendrank = iditr - readdispls.cbegin();
        size_t reqstart = std::max(static_cast<size_t>(*iditr++), startid);
        size_t reqend = std::min(static_cast<size_t>(*iditr), startid + count);
        if (sendrank != myrank) mynbrs.push_back({reqstart, reqend-reqstart, sendrank, myrank, rc_flag});
    }
}

std::ostream& operator<<(std::ostream& stream, const NbrData& o)
{
    stream << "rank " << o.destrank+1 << " requests " << "[" << o.startid << ".." << o.startid+o.numseqs-1 << "] from rank " << o.sendrank+1 << " (rc=" << o.rcflag << ")";
    return stream;
}

void DistributedFastaData::allgather_neighbors()
{
    Grid commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    /*
     * TODO: update this comment
     *
     * allnbrs stores all triples (NbrData is an (int,int,int)) (sendrank, destrank, rcflag)
     * where the sendrank is the rank of the processor sending its
     * FastaData, the destrank is the rank of processor receiving
     * some FastaData, and rcflag is 0 or 1 depending on whether the
     * data being sent is from the row slice or the column slice.
     *
     * mynbrs are the triples I am receiving. We will allgather
     * them so that everyone knows everything eventually.
     */
    std::vector<NbrData> mynbrs, allnbrs;

    /*
     * Add all the row data I need to receive to mynbrs.
     */
    findnbrs(mynbrs, rowstartid, numrowreads, 0);

    /*
     * Add all the column data I need to receive to mynbrs (if not a diagonal processor).
     */
    if (myrowid != mycolid) findnbrs(mynbrs, colstartid, numcolreads, 1);

    /*
     * Allgather NbrData.
     */
    MPI_Count_type allnbrcount, mynbrcount = mynbrs.size();
    std::vector<MPI_Count_type> nbrcounts(nprocs);
    std::vector<MPI_Displ_type> nbrdispls(nprocs);

    nbrcounts[myrank] = mynbrcount;
    MPI_ALLGATHER(MPI_IN_PLACE, 1, MPI_COUNT_TYPE, nbrcounts.data(), 1, MPI_COUNT_TYPE, comm);

    std::exclusive_scan(nbrcounts.begin(), nbrcounts.end(), nbrdispls.begin(), static_cast<MPI_Displ_type>(0));
    allnbrcount = nbrcounts.back() + nbrdispls.back();

    MPI_Datatype nbrdtype;
    int blklens[5] = {1,1,1,1,1};
    MPI_Aint displs[5] = {offsetof(NbrData, startid),
                          offsetof(NbrData, numseqs),
                          offsetof(NbrData, sendrank),
                          offsetof(NbrData, destrank),
                          offsetof(NbrData, rcflag)};

    MPI_Datatype types[5] = {MPI_SIZE_T, MPI_SIZE_T, MPI_INT, MPI_INT, MPI_UNSIGNED_SHORT};
    MPI_Type_create_struct(5, blklens, displs, types, &nbrdtype);
    MPI_Type_commit(&nbrdtype);

    allnbrs.resize(allnbrcount);
    MPI_ALLGATHERV(mynbrs.data(), mynbrcount, nbrdtype, allnbrs.data(), nbrcounts.data(), nbrdispls.data(), nbrdtype, comm);

    MPI_Type_free(&nbrdtype);
}

void DistributedFastaData::exchange_reads()
{
    Grid commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();
}
