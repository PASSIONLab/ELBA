#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"
#include "DnaBuffer.hpp"

class FastaIndex
{
public:
    typedef struct { size_t len, pos, bases; } Record;

    FastaIndex(const std::string& fasta_fname, Grid commgrid);

    Grid getcommgrid() const { return commgrid; }
    std::string get_fasta_fname() const { return fasta_fname; }
    std::string get_faidx_fname() const { return fasta_fname + ".fai"; }

    size_t gettotrecords() const { return readdispls.back(); }
    size_t getreadcount(size_t i) const { return static_cast<size_t>(readcounts[i]); }
    size_t getreaddispl(size_t i) const { return static_cast<size_t>(readdispls[i]); }
    size_t getmyreadcount() const { return getreadcount(commgrid->GetRank()); }
    size_t getmyreaddispl() const { return getreaddispl(commgrid->GetRank()); }
    size_t totbases() const { return std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return sum + record.len; }); }
    size_t maxlen() const { return std::accumulate(myrecords.begin(), myrecords.end(), static_cast<size_t>(0), [](size_t sum, const auto& record) { return std::max(sum, record.len); }); }

    const std::vector<Record>& getmyrecords() const { return myrecords; }
    const std::vector<MPI_Count_type> getreadcounts() const { return readcounts; }
    const std::vector<MPI_Displ_type> getreaddispls() const { return readdispls; }

    std::shared_ptr<DnaBuffer> getmydna() const;

    static Record get_faidx_record(const std::string& line);

private:
    Grid commgrid;
    std::vector<Record> myrecords; /* records for the reads local processor is responsible for */
    std::vector<Record> rootrecords; /* records for all the reads across all processors (parsed by root rank in constructor) */
    std::vector<MPI_Count_type> readcounts; /* number of reads assigned to each processor. Each processor gets a copy. |readcounts| == nprocs */
    std::vector<MPI_Displ_type> readdispls; /* displacement counts for reads across all processors. Each processor gets a copy. |readdispls| == nprocs+1 */
    std::string fasta_fname; /* FASTA file name */

    void getpartition(std::vector<MPI_Count_type>& sendcounts);
};

#endif
