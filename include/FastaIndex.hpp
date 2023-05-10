#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"
#include "DnaBuffer.hpp"

class FastaIndex
{
public:
    typedef struct { size_t len, pos, bases; } Record;

    FastaIndex(const std::string& fasta_fname, std::shared_ptr<CommGrid> commgrid);

    std::shared_ptr<CommGrid> getcommgrid() const { return commgrid; }
    std::string get_fasta_fname() const { return fasta_fname; }
    std::string get_faidx_fname() const { return fasta_fname + ".fai"; }

    size_t gettotrecords() const { return readdispls.back(); }
    size_t getreadcount(size_t i) const { return static_cast<size_t>(readcounts[i]); }
    size_t getreaddispl(size_t i) const { return static_cast<size_t>(readdispls[i]); }
    size_t getmyreadcount() const { return getreadcount(commgrid->GetRank()); }
    size_t getmyreaddispl() const { return getreaddispl(commgrid->GetRank()); }

    std::vector<size_t> getmyreadlens() const;
    std::vector<size_t> getallreadlens() const;

    const std::vector<Record>& getmyrecords() const { return myrecords; }
    const std::vector<MPI_Count_type> getreadcounts() const { return readcounts; }
    const std::vector<MPI_Displ_type> getreaddispls() const { return readdispls; }

    DnaBuffer getmydna() const;
    void log(const DnaBuffer& buffer) const;

    static Record get_faidx_record(const std::string& line, std::string& name);

    std::vector<std::string> bcastnames();

private:
    std::shared_ptr<CommGrid> commgrid;
    std::vector<Record> myrecords; /* records for the reads local processor is responsible for */
    std::vector<Record> rootrecords; /* records for all the reads across all processors (parsed by root rank in constructor) */
    std::vector<MPI_Count_type> readcounts; /* number of reads assigned to each processor. Each processor gets a copy. |readcounts| == nprocs */
    std::vector<MPI_Displ_type> readdispls; /* displacement counts for reads across all processors. Each processor gets a copy. |readdispls| == nprocs+1 */
    std::string fasta_fname; /* FASTA file name */

    std::vector<std::string> rootnames;

    void getpartition(std::vector<MPI_Count_type>& sendcounts);
};

#endif
