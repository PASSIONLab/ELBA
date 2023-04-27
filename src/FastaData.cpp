#include "FastaData.hpp"
#include "DnaSeq.hpp"
#include "Logger.hpp"

FastaData::FastaData(FIndex index) : index(index), idxtag(index->getcommgrid()->GetRank())
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

    readlens.reserve(numreads);
    byteoffsets.reserve(numreads+1);

    size_t totbases = 0, maxlen = 0;
    for (auto itr = records.cbegin(); itr != records.cend(); ++itr)
    {
        totbases += itr->len;
        maxlen = std::max(itr->len, maxlen);
        readlens.push_back(itr->len);
    }

    size_t bufbound = FastaData::computebufbound(totbases, numreads);
    buf = new uint8_t[bufbound];

    char *tmpbuf = new char[maxlen];
    size_t byteoffset = 0;

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

        size_t numbytes = FastaData::bytesneeded(itr->len);
        size_t overhang = (numbytes * 4) - itr->len;
        assert(overhang < 4);

        size_t b = 0;
        char const *sb = tmpbuf;
        uint8_t *bytes = buf + byteoffset;
        uint8_t byte;

        while (b < numbytes)
        {
            byte = 0;
            int left = (b != numbytes-1? 4 : 4-overhang);

            for (int i = 0; i < left; ++i)
            {
                uint8_t code = DnaSeq::getcharcode(sb[i]);
                uint8_t shift = code << (6 - (2*i));
                byte |= shift;
            }

            bytes[b++] = byte;
            sb += 4;
        }

        byteoffsets.push_back(byteoffset);
        byteoffset += b;
    }

    byteoffsets.push_back(byteoffset); /* convenient */
    assert(bufbound >= byteoffset); /* one computed from principals, the other computed via above algorithm, so if this fails then algorithm above is wrong */

    elapsed += MPI_Wtime();
    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.Flush("FASTA parsing rates (FastaData):");

    delete[] tmpbuf;
    delete[] readbuf;
}

std::string FastaData::getsequence(size_t localid) const
{
    /*
     * Probably never use this function, however it provides
     * an understanding of how the data is actually stored.
     */

    size_t len = readlens[localid];
    uint8_t const *bytes = buf + byteoffsets[localid];
    std::vector<char> s(len);

    for (size_t i = 0; i < len; ++i)
    {
        int code = (*bytes >> (6 - (2 * (i%4)))) & 3;
        s[i] = DnaSeq::getcodechar(code);

        if ((i+1) % 4 == 0)
            bytes++;
    }

    return std::string(s.begin(), s.end());
}

void FastaData::log() const
{
    Logger logger(index->getcommgrid());
    assert(index->getnumrecords() == readlens.size());
    size_t numreads = index->getnumrecords();
    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0));
    double avglen = static_cast<double>(totbases) / numreads;
    size_t firstid = getfirstid();
    logger() << "idxtag " << idxtag << " stores reads " << Logger::readrangestr(firstid, numreads) << ". ~" << std::fixed << std::setprecision(2) << avglen << " nts/read. (" << static_cast<double>(getbufsize()) / (1024.0 * 1024.0) << " Mbs compressed) == (" << getbufsize() << " bytes)";
    logger.Flush("FASTA sequence distribution (FastaData):");
}

size_t FastaData::computebufbound(size_t totbases, size_t numreads)
{
    /*
     * For a sequence of length l, floor((l+3)/4) bytes are needed to encode it. Let
     * L = l[0] + l[1] + ... + l[N-1], where l[i] is the length of the ith sequence, N is
     * the total number of sequences, and hence L is the total sum of all the sequence
     * lengths. Then the total number of bytes needed for the write buffer is
     *
     * Sum{0 <= i <= N-1}[floor((l[i]+3)/4)] <= (1/4) * Sum{0 <= i <= N-1}[l[i] + 4]
     *                                       <= (1/4) * (L + 4N)
     *                                        = (L/4) + N
     */
    return static_cast<size_t>(std::ceil((totbases/4.0) + numreads));
}

size_t FastaData::bytesneeded(size_t numbases)
{
    return (numbases + 3) / 4;
}
