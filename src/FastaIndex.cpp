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

Record FastaIndex::get_faidx_record(const std::string& line, std::string& name)
{
    /*
     * Read a line from a FASTA index file into a record object.
     */
    Record record;
    std::istringstream(line) >> name >> record.len >> record.pos >> record.bases;
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

    /*
     * Coming up with the optimal partitioning of sequences weighted by their length
     * is NP-hard, and I can't imagine that it is NP-complete also (how would you verify
     * in polynomial time?). Tried looking for the name of this optimization problem
     * but couldn't find (didn't look too hard because its not too important for this).
     * It's basically a variation on multiway number partition, except where the divisions
     * must be ordered. The following seems like a fine approxmimation. Note that the
     * the last processor tends to get more than the average amount of data.
     */

    size_t readid = 0;

    for (int i = 0; i < nprocs-1; ++i)
    {
        size_t basessofar = 0;
        size_t startid = readid;

        /*
         * Keep going through reads until the next one puts us over
         * the average bases per processor.
         */
        while (readid < numreads && basessofar + rootrecords[readid].len < avgbasesperproc)
        {
            basessofar += rootrecords[readid].len;
            readid++;
        }

        size_t readssofar = readid - startid;
        assert(readssofar >= 1); /* TODO: come up with a recovery strategy here */

        sendcounts[i] = readssofar;
    }

    /*
     * Last processor gets the remaining reads.
     */
    sendcounts.back() = numreads - readid;
}

FastaIndex::FastaIndex(const std::string& fasta_fname, std::shared_ptr<CommGrid> commgrid) : commgrid(commgrid), fasta_fname(fasta_fname)
{
    int nprocs = commgrid->GetSize();
    int myrank = commgrid->GetRank();
    MPI_Comm comm = commgrid->GetWorld();
    readcounts.resize(nprocs);

    /*
     * Root processor responsible for reading and parsing FASTA
     * index file "{fasta_fname}.fai" into one record per sequence.
     */
    if (myrank == 0)
    {
        std::string line, name;
        std::ifstream filestream(get_faidx_fname());

        while (std::getline(filestream, line))
        {
            rootrecords.push_back(get_faidx_record(line, name));
            rootnames.push_back(name);
        }

        filestream.close();

        /*
         * Compute load-balanced read partitioning on the root processor.
         */
        getpartition(readcounts);
    }

    /*
     * All processors get a copy of the read counts.
     */
    MPI_BCAST(readcounts.data(), nprocs, MPI_COUNT_TYPE, 0, comm);

    /*
     * And the read displacements.
     */
    readdispls.resize(nprocs);
    std::exclusive_scan(readcounts.begin(), readcounts.end(), readdispls.begin(), static_cast<MPI_Displ_type>(0));

    /*
     * It is useful for the displacements to store the total number
     * of reads in the last position.
     */
    readdispls.push_back(readdispls.back() + readcounts.back());

    /*
     * To prevent confusion for the reader: the broadcasting of readcounts
     * to every processor and the parallel "re"-computation of readdispls
     * on each processor is not necessary for performing the scatter operation
     * below, however we still do it because every processor will need to know
     * those things later.
     */
    myrecords.resize(readcounts[myrank]);

    /*
     * Each record is represented with three numbers (read length,
     * FASTA position, FASTA line width), so we create an MPI datatype
     * to communicate each record as a single unit.
     */
    MPI_Datatype faidx_dtype_t;
    MPI_Type_contiguous(3, MPI_SIZE_T, &faidx_dtype_t);
    MPI_Type_commit(&faidx_dtype_t);

    /*
     * Scatter the records according to the load-balanced read partitioning.
     */
    MPI_SCATTERV(rootrecords.data(), readcounts.data(), readdispls.data(), faidx_dtype_t, myrecords.data(), readcounts[myrank], faidx_dtype_t, 0, comm);

    MPI_Type_free(&faidx_dtype_t);

    #if LOG_LEVEL >= 2
    Logger logger(commgrid);
    size_t mytotbases = std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; });
    size_t totbases;
    MPI_ALLREDUCE(&mytotbases, &totbases, 1, MPI_SIZE_T, MPI_SUM, comm);
    double percent_proportion = (static_cast<double>(mytotbases) / totbases) * 100.0;
    logger() << " is responsible for sequences " << Logger::readrangestr(readdispls[myrank], readcounts[myrank]) << " (" << mytotbases << " nucleotides, " << std::fixed << std::setprecision(3) << percent_proportion << "%)";
    logger.Flush("Fasta index construction:");
    #endif
}

std::vector<size_t> FastaIndex::getmyreadlens() const
{
    /*
     *  Because we store sequence information using "array of structs"
     *  instead of "structs of arrays", if we want a linear array
     *  of read lengths we have to unpack them.
     */
    std::vector<size_t> readlens(getmyreadcount());
    std::transform(myrecords.cbegin(), myrecords.cend(), readlens.begin(), [](const auto& record) { return record.len; });
    return readlens;
}

DnaBuffer FastaIndex::getmydna() const
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    /*
     * Allocate local sequence buffer.
     */
    auto readlens = getmyreadlens(); /* vector of local read lengths */
    size_t bufsize = DnaBuffer::computebufsize(readlens); /* minimum number of bytes needed to 2-bit encode all the local reads */
    DnaBuffer dnabuf(bufsize); /* initialize dnabuf by allocating @bufsize bytes */
    size_t numreads = readlens.size(); /* number of local reads */

    MPI_Offset startpos; /* the FASTA position that starts my local chunk of reads */
    MPI_Offset endpos; /* the FASTA position that ends my local chunk of reads (exclusive) */
    MPI_Offset filesize; /* the total size of the FASTA */
    MPI_Offset readbufsize; /* endpos - startpos */
    MPI_File fh;

    /*
     * FASTA file will be read using MPI collective I/O.
     */
    MPI_File_open(comm, get_fasta_fname().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    MPI_File_get_size(fh, &filesize);

    /*
     * Get start and end coordinates within FASTA of the sequences
     * this processor requires.
     */
    startpos = myrecords.front().pos;
    endpos = myrecords.back().pos + myrecords.back().len + (myrecords.back().len / myrecords.back().bases);
    if (endpos > filesize) endpos = filesize;

    /*
     * Allocate a char buffer to read my FASTA chunk into to,
     * and then do the reading.
     */
    readbufsize = endpos - startpos;
    std::unique_ptr<char[]> readbuf(new char[readbufsize]);

    /*
     * Every processor that calls @getmydna reads in its assigned chunk of data
     * into its temporary local storage buffer (meant for raw contents of file).
     */
    MPI_FILE_READ_AT_ALL(fh, startpos, &readbuf[0], readbufsize, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);

    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0), std::plus<size_t>{});
    size_t maxlen = *std::max_element(readlens.begin(), readlens.end());

    /*
     * ASCII sequences are first read into this temporary char buffer
     * before they are compressed into the sequence buffer.
     */
    std::unique_ptr<char[]> tmpbuf(new char[maxlen]);

    MPI_Barrier(comm);
    double elapsed = -MPI_Wtime();

    /*
     * Go through each local FASTA record.
     */
    for (auto itr = myrecords.cbegin(); itr != myrecords.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *writeptr = &tmpbuf[0];

        /*
         * Read ASCII FASTA sequence into the temoprary buffer.
         */
        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(writeptr, &readbuf[chunkpos + locpos], cnt);
            writeptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        /*
         * DnaBuffer automatically 2-bit encodes the ASCII sequence
         * and pushes it onto its local stack of sequences.
         */
        dnabuf.push_back(&tmpbuf[0], itr->len);
    }

    elapsed += MPI_Wtime();

    #if LOG_LEVEL >= 2
    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger(commgrid);
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.Flush("FASTA parsing rates (DnaBuffer):");
    #endif

    return dnabuf;
}

std::vector<std::string> FastaIndex::bcastnames()
{
    /*
     * In order to print out PAF files, we need the names for each
     * read sequence. It's easiest to just broadcast all the
     * read names to every process (although it's not the most efficient).
     */

    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    MPI_Count_type numchars; /* total number of characters needed for every read name */
    MPI_Count_type numreads = gettotrecords(); /* total number of reads (and hence names) */
    std::vector<char> flatbuf; /* this will store a flattened buffer of all the characters linearly ordered from all the names */
    std::vector<size_t> namelens(numreads); /* vector of read name lengths for unflattening the flatbuf */

    /*
     * At the start, the root process stores all the read names and the other
     * processes don't know any of the names. So the root process computes the
     * the total number of characters needed and then broadcasts them.
     */
    if (myrank == 0)
    {
        numchars = std::accumulate(rootnames.begin(), rootnames.end(), static_cast<MPI_Count_type>(0), [](MPI_Count_type sum, const std::string& s) { return sum + s.size(); });
    }

    MPI_BCAST(&numchars, 1, MPI_COUNT_TYPE, 0, comm);

    /*
     * The root process goes through each name in the rootnames vector
     * and appends them onto the flat buffer. For example, if
     * rootnames = {'read1', read2', 'read3'}, then flatbuf = "read1read2read3".
     */
    if (myrank == 0)
    {
        flatbuf.reserve(numchars);
        auto itr = rootnames.begin();

        for (size_t i = 0; i < static_cast<size_t>(numreads); ++i)
        {
            flatbuf.insert(flatbuf.end(), itr->begin(), itr->end()); /* appends name on the end of flatbuf */
            namelens[i] = itr->size(); /* record the length of the name */
            ++itr;
        }
    }
    else
    {
        /*
         * Non-root processes use the broadcasted flatbuf length to allocate
         * the memory needed to store its copy of flatbuf.
         */
        flatbuf.resize(numchars);
    }

    assert(flatbuf.size() == numchars);

    /*
     * Broadcast the name lengths and the flat buffer.
     */
    MPI_BCAST(namelens.data(), numreads, MPI_SIZE_T, 0, comm);
    MPI_BCAST(flatbuf.data(), numchars, MPI_CHAR, 0, comm);

    /*
     * Root process can return immediately because its
     * vector of names has already been computed.
     */
    if (myrank == 0)
    {
        return rootnames;
    }

    /*
     * Non-root processes have to unpack flatbuf using the
     * name lengths before returning.
     */
    std::vector<std::string> recvnames;
    recvnames.reserve(numreads);

    auto itr = flatbuf.begin();

    for (size_t i = 0; i < static_cast<size_t>(numreads); ++i)
    {
        recvnames.emplace_back(itr, itr + namelens[i]);
        itr += namelens[i];
    }

    return recvnames;
}

void FastaIndex::log(const DnaBuffer& buffer) const
{
    Logger logger(commgrid);
    assert(getmyreadcount() == buffer.size());
    auto readlens = getmyreadlens();
    size_t numreads = getmyreadcount();
    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0), std::plus<size_t>{});
    double avglen = static_cast<double>(totbases) / numreads;
    size_t firstid = getmyreaddispl();
    logger() << " stores " << Logger::readrangestr(firstid, numreads) << ". ~" << std::fixed << std::setprecision(2) << avglen << " nts/read. (" << static_cast<double>(buffer.getbufsize()) / (1024.0 * 1024.0) << " Mbs compressed) == (" << buffer.getbufsize() << " bytes)";
    logger.Flush("FASTA sequence storage (DnaBuffer):");
}
