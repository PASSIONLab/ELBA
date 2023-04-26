#include "FastaIndex.hpp"
#include "Logger.hpp"
#include <cstring>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>

typedef typename FastaIndex::faidx_record_t faidx_record_t;

faidx_record_t GetFaidxRecord(const std::string& line, std::vector<std::string>& names)
{
    std::string name;
    faidx_record_t record;
    std::istringstream(line) >> name >> record.len >> record.pos >> record.bases;
    names.push_back(name);
    return record;
}

MPI_Count_type FastaIndex::get_idbalanced_partition(std::vector<MPI_Count_type>& sendcounts)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    size_t numreads = GetNumRecords();
    sendcounts.resize(nprocs);

    MPI_Count_type readsperproc = numreads / nprocs;

    std::fill_n(sendcounts.begin(), nprocs-1, readsperproc);

    sendcounts.back() = numreads - (nprocs-1) * readsperproc;

    return numreads;
}

MPI_Count_type FastaIndex::get_membalanced_partition(std::vector<MPI_Count_type>& sendcounts)
{
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    size_t numreads = GetNumRecords();
    sendcounts.resize(nprocs);

    size_t totbases = std::accumulate(records.begin(), records.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });
    double avgbasesperproc = static_cast<double>(totbases) / nprocs;

    size_t readid = 0;

    for (int i = 0; i < nprocs-1; ++i)
    {
        size_t basessofar = 0;
        size_t startid = readid;

        while (readid < numreads && basessofar + records[readid].len < avgbasesperproc)
        {
            basessofar += records[readid].len;
            readid++;
        }

        size_t readssofar = readid - startid;
        assert(readssofar >= 1);

        sendcounts[i] = readssofar;
    }

    sendcounts.back() = numreads - readid;
    return numreads;
}

FastaIndex::FastaIndex(const std::string& fasta_fname, Grid commgrid, bool membalanced) : commgrid(commgrid), fasta_fname(fasta_fname)
{
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();

    std::vector<MPI_Count_type> sendcounts; /* MPI_Scatterv sendcounts for faidx_record_t records (root only) */
    std::vector<MPI_Displ_type> displs;     /* MPI_Scatterv displs for faidx_record_t records (root only)     */
    MPI_Count_type recvcount;               /* MPI_Scatterv recvcount for faidx_record_t records              */

    std::vector<std::string> root_names;

    if (myrank == 0)
    {
        std::string line;
        std::ifstream filestream(GetFaidxFilename());

        while (std::getline(filestream, line))
        {
            records.push_back(GetFaidxRecord(line, root_names));
        }

        filestream.close();

        MPI_Count_type num_records;

        num_records = membalanced? get_membalanced_partition(sendcounts) : get_idbalanced_partition(sendcounts);

        displs.resize(nprocs);
        std::exclusive_scan(sendcounts.begin(), sendcounts.end(), displs.begin(), static_cast<MPI_Displ_type>(0));
    }

    /*
     * Root process tells each process how many faidx_record_t records it will be sent.
     */
    MPI_SCATTER(sendcounts.data(), 1, MPI_COUNT_TYPE, &recvcount, 1, MPI_COUNT_TYPE, 0, commgrid->GetWorld());

    myrecords.resize(recvcount);

    MPI_Datatype faidx_dtype_t;
    MPI_Type_contiguous(3, MPI_SIZE_T, &faidx_dtype_t);
    MPI_Type_commit(&faidx_dtype_t);
    MPI_SCATTERV(records.data(), sendcounts.data(), displs.data(), faidx_dtype_t, myrecords.data(), recvcount, faidx_dtype_t, 0, commgrid->GetWorld());
    MPI_Type_free(&faidx_dtype_t);

    size_t mynumreads = myrecords.size();
    size_t mytotbases = std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });

    Logger logger(commgrid);
    logger() << " responsible for " << mynumreads << " reads with a total of " << mytotbases << " nucleotides according to " << GetFaidxFilename();
    logger.Flush("FastaIndex::FastaIndex");
}

