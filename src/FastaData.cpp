#include "FastaData.hpp"
#include "DnaSeq.hpp"
#include "Logger.hpp"

FastaData::FastaData(FIndex index)
{
    /*
     * Will skip convention of adding 'my' to the front of variables referring
     * to local processor information since that is implied by the nature
     * of this class' purpose.
     */
    size_t numreads;
    char *readbuf;
    MPI_Offset startpos, endpos, filesize, readbufsize;
    MPI_File fh;

    Grid commgrid = index->getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    const auto& records = index->getmyrecords();
    numreads = records.size();
    startpos = records.front().pos;
    endpos = records.back().pos + records.back().len + (records.back().len / records.back().bases);
    MPI_File_open(comm, index->get_fasta_fname().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &filesize);
    if (endpos > filesize) endpos = filesize;
    readbufsize = endpos - startpos;
    readbuf = new char[readbufsize];
    MPI_FILE_READ_AT_ALL(fh, startpos, readbuf, readbufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    size_t totbases = index->totbases();
    size_t maxlen = index->maxlen();

    buffer.reset(new DnaBuffer(totbases, numreads));

    char *tmpbuf = new char[maxlen];

    MPI_Barrier(comm);
    double elapsed = -MPI_Wtime();

    for (auto itr = records.cbegin(); itr != records.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *writeptr = tmpbuf;

        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(writeptr, &readbuf[chunkpos + locpos], cnt);
            writeptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        /* ASCII sequence is now stored in tmpbuf */

        sequences.emplace_back(tmpbuf, itr->len, *buffer);
    }

    elapsed += MPI_Wtime();
    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.Flush("FASTA parsing rates (FastaData):");

    delete[] tmpbuf;
    delete[] readbuf;
}


void FastaData::log(FIndex index) const
{
    Logger logger(index->getcommgrid());
    assert(index->getmyreadcount() == buffer->getnumseqs());
    size_t numreads = index->getmyreadcount();
    size_t totbases = index->totbases();
    double avglen = static_cast<double>(totbases) / numreads;
    size_t firstid = index->getmyreaddispl();
    logger() << " stores " << Logger::readrangestr(firstid, numreads) << ". ~" << std::fixed << std::setprecision(2) << avglen << " nts/read. (" << static_cast<double>(buffer->getbufsize()) / (1024.0 * 1024.0) << " Mbs compressed) == (" << buffer->getbufsize() << " bytes)";
    logger.Flush("FASTA sequence distribution (FastaData):");
}
