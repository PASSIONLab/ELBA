#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "common.h"

class FastaIndex
{
public:
    typedef struct faidx_record { size_t len, pos, bases; } faidx_record_t;

    FastaIndex(const std::string& fasta_fname, Grid commgrid, bool membalanced = false);

    Grid getcommgrid() const { return commgrid; }

    size_t GetNumRecords() const { return records.size(); }

    const std::vector<faidx_record_t>& getmyrecords() const { return myrecords; }

    std::string GetFastaFilename() const { return fasta_fname; }
    std::string GetFaidxFilename() const { return fasta_fname + ".fai"; }

private:
    Grid commgrid;
    std::vector<faidx_record_t> myrecords, records;
    std::vector<std::string> names;
    std::string fasta_fname;

    MPI_Count_type get_idbalanced_partition(std::vector<MPI_Count_type>& sendcounts);
    MPI_Count_type get_membalanced_partition(std::vector<MPI_Count_type>& sendcounts);
};

#endif
