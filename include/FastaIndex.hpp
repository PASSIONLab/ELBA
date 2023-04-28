#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"

class FastaIndex
{
public:
    typedef struct { size_t len, pos, bases; } Record;

    FastaIndex(const std::string& fasta_fname, Grid commgrid);

    Grid getcommgrid() const { return commgrid; }
    size_t getnumrecords() const { return readcounts[commgrid->GetRank()]; }
    size_t gettotrecords() const { return readdispls.back(); }
    std::string get_fasta_fname() const { return fasta_fname; }
    std::string get_faidx_fname() const { return fasta_fname + ".fai"; }
    size_t getsomefirstid(int rank) const { return static_cast<size_t>(readdispls[rank]); }
    size_t getmyfirstid() const { return getsomefirstid(commgrid->GetRank()); }
    size_t getreadcount(int rank) const { return static_cast<size_t>(readcounts[rank]); }
    std::vector<int> collectowners(size_t startid, size_t count) const;
    static Record get_faidx_record(const std::string& line);

    const std::vector<Record>& getmyrecords() const { return myrecords; }
    const std::vector<MPI_Displ_type>& getreaddispls() const { return readdispls; }

private:
    Grid commgrid;
    std::vector<Record> myrecords; /* records for the reads local processor is responsible for */
    std::vector<Record> rootrecords; /* records for all the reads across all processors (parsed by root rank in constructor) */
    std::vector<MPI_Count_type> readcounts; /* number of reads assigned to each processor. Each processor gets a copy. |readcounts| == nprocs */
    std::vector<MPI_Displ_type> readdispls; /* displacement counts for reads across all processors. Each processor gets a copy. |readdispls| == nprocs+1 */
    std::string fasta_fname; /* FASTA file name */

    void get_idbalanced_partition(std::vector<MPI_Count_type>& sendcounts);
    void get_membalanced_partition(std::vector<MPI_Count_type>& sendcounts);
};

using FIndex = std::shared_ptr<FastaIndex>;

#endif
