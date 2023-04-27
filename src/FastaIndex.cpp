#include "FastaIndex.hpp"
#include "Logger.hpp"
#include <cstring>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>

using Record = typename FastaIndex::Record;

Record FastaIndex::get_faidx_record(const std::string& line)
{
    std::string dummy;
    Record record;
    std::istringstream(line) >> dummy >> record.len >> record.pos >> record.bases;
    return record;
}
void FastaIndex::get_idbalanced_partition(std::vector<MPI_Count_type>& sendcounts)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    size_t numreads = rootrecords.size();
    sendcounts.resize(nprocs);

    MPI_Count_type readsperproc = numreads / nprocs;

    std::fill_n(sendcounts.begin(), nprocs-1, readsperproc);

    sendcounts.back() = numreads - (nprocs-1) * readsperproc;
}

void FastaIndex::get_membalanced_partition(std::vector<MPI_Count_type>& sendcounts)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    assert(sendcounts.size() == nprocs);

    size_t totbases = std::accumulate(rootrecords.begin(), rootrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });
    size_t numreads = rootrecords.size();
    double avgbasesperproc = static_cast<double>(totbases) / nprocs;

    size_t readid = 0;

    for (int i = 0; i < nprocs-1; ++i)
    {
        size_t basessofar = 0;
        size_t startid = readid;

        while (readid < numreads && basessofar + rootrecords[readid].len < avgbasesperproc)
        {
            basessofar += rootrecords[readid].len;
            readid++;
        }

        size_t readssofar = readid - startid;
        assert(readssofar >= 1);

        sendcounts[i] = readssofar;
    }

    sendcounts.back() = numreads - readid;
}

FastaIndex::FastaIndex(const std::string& fasta_fname, Grid commgrid) : commgrid(commgrid), fasta_fname(fasta_fname)
{
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    readcounts.resize(nprocs);

    if (myrank == 0)
    {
        std::string line;
        std::ifstream filestream(get_faidx_fname());
        while (std::getline(filestream, line)) rootrecords.push_back(get_faidx_record(line));
        filestream.close();
        get_membalanced_partition(readcounts);
    }

    MPI_BCAST(readcounts.data(), nprocs, MPI_COUNT_TYPE, 0, comm);

    readdispls.resize(nprocs);
    std::exclusive_scan(readcounts.begin(), readcounts.end(), readdispls.begin(), static_cast<MPI_Displ_type>(0));
    readdispls.push_back(readdispls.back() + readcounts.back());
    assert(readdispls[nprocs] - readdispls[nprocs-1] == readcounts.back());

    /*
     * To prevent confusion for the reader: the broadcasting of readcounts
     * to every processor and the the parallel "re"-computation of readdispls
     * on each processor is not necessary for performing the scatter operation
     * below, however we still do it because every processor will need to know
     * those things later.
     */

    myrecords.resize(readcounts[myrank]);

    MPI_Datatype faidx_dtype_t;
    MPI_Type_contiguous(3, MPI_SIZE_T, &faidx_dtype_t);
    MPI_Type_commit(&faidx_dtype_t);
    MPI_SCATTERV(rootrecords.data(), readcounts.data(), readdispls.data(), faidx_dtype_t, myrecords.data(), readcounts[myrank], faidx_dtype_t, 0, comm);
    MPI_Type_free(&faidx_dtype_t);

    Logger logger(commgrid);
    size_t mytotbases = std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });
    size_t totbases;
    MPI_ALLREDUCE(&mytotbases, &totbases, 1, MPI_SIZE_T, MPI_SUM, comm);
    double percent_proportion = (static_cast<double>(mytotbases) / totbases) * 100.0;
    logger() << " is responsible for sequences " << Logger::readrangestr(readdispls[myrank], readcounts[myrank]) << " (" << mytotbases << " nucleotides, " << std::fixed << std::setprecision(3) << percent_proportion << "%)";
    logger.Flush("Fasta index construction:");
}
