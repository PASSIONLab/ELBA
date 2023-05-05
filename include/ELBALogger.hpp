#ifndef ELBA_LOGGER_H_
#define ELBA_LOGGER_H_

#include "common.h"
#include "KmerOps.hpp"
#include "SharedSeeds.hpp"
#include "Overlap.hpp"

struct ELBALogger
{
    std::string prefix;
    MPI_Comm comm;
    bool isroot;

    ELBALogger(const std::string& prefix, MPI_Comm comm) : prefix(prefix), comm(comm)
    {
        int myrank;
        MPI_Comm_rank(comm, &myrank);
        isroot = (myrank == 0);
    }

    void log_kmer_matrix(const CT<PosInRead>::PSpParMat& A);
    void log_seed_matrix(const CT<SharedSeeds>::PSpParMat& B);
    void log_overlap_matrix(const CT<Overlap>::PSpParMat& R);

    std::string getmatfname(const std::string matname);
};

#endif
