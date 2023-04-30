#include "FastaIndex.hpp"
#include "Logger.hpp"
#include <cstring>
#include <iterator>
#include <algorithm>
#include <functional>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <memory>

using Record = typename FastaIndex::Record;

Record FastaIndex::get_faidx_record(const std::string& line)
{
    std::string dummy;
    Record record;
    std::istringstream(line) >> dummy >> record.len >> record.pos >> record.bases;
    return record;
}

void FastaIndex::getpartition(std::vector<MPI_Count_type>& sendcounts)
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
        getpartition(readcounts);
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

std::vector<size_t> FastaIndex::getmyreadlens() const
{
    std::vector<size_t> readlens(getmyreadcount());
    std::transform(myrecords.cbegin(), myrecords.cend(), readlens.begin(), [](const auto& record) { return record.len; });
    return readlens;
}

std::shared_ptr<DnaBuffer> FastaIndex::getmydna() const
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    auto readlens = getmyreadlens();
    size_t bufsize = DnaBuffer::computebufsize(readlens);
    auto dnabuf = std::make_shared<DnaBuffer>(bufsize);

    size_t numreads;
    MPI_Offset startpos, endpos, filesize, readbufsize;
    MPI_File fh;

    numreads = myrecords.size();
    startpos = myrecords.front().pos;
    endpos = myrecords.back().pos + myrecords.back().len + (myrecords.back().len / myrecords.back().bases);

    MPI_File_open(comm, get_fasta_fname().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &filesize);
    if (endpos > filesize) endpos = filesize;
    readbufsize = endpos - startpos;
    std::unique_ptr<char[]> readbuf(new char[readbufsize]);
    MPI_FILE_READ_AT_ALL(fh, startpos, &readbuf[0], readbufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0), std::plus<size_t>{});
    size_t maxlen = *std::max_element(readlens.begin(), readlens.end());

    std::unique_ptr<char[]> tmpbuf(new char[maxlen]);

    MPI_Barrier(comm);
    double elapsed = -MPI_Wtime();

    for (auto itr = myrecords.cbegin(); itr != myrecords.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *writeptr = &tmpbuf[0];

        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(writeptr, &readbuf[chunkpos + locpos], cnt);
            writeptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        dnabuf->push_back(&tmpbuf[0], itr->len);
    }

    elapsed += MPI_Wtime();
    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.Flush("FASTA parsing rates (FastaData):");

    return dnabuf;
}

void FastaIndex::log(std::shared_ptr<DnaBuffer> buffer) const
{
    Logger logger(commgrid);
    assert(getmyreadcount() == buffer->size());
    auto readlens = getmyreadlens();
    size_t numreads = getmyreadcount();
    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0), std::plus<size_t>{});
    double avglen = static_cast<double>(totbases) / numreads;
    size_t firstid = getmyreaddispl();
    logger() << " stores " << Logger::readrangestr(firstid, numreads) << ". ~" << std::fixed << std::setprecision(2) << avglen << " nts/read. (" << static_cast<double>(buffer->getbufsize()) / (1024.0 * 1024.0) << " Mbs compressed) == (" << buffer->getbufsize() << " bytes)";
    logger.Flush("FASTA sequence storage (DnaBuffer):");
}
