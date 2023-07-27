#include "ELBALogger.hpp"
#include <iostream>
#include <sstream>

void ELBALogger::log_kmer_matrix(CT<PosInRead>::PSpParMat& A)
{
    #if LOG_LEVEL >= 1
    size_t numreads = A.getnrow();
    size_t numkmers = A.getncol();
    size_t numseeds = A.getnnz();

    if (isroot) std::cout << "K-mer matrix A has " << numreads << " rows (readids), " << numkmers << " columns (k-mers), and " << numseeds << " nonzeros (k-mer seeds)\n" << std::endl;

    #if LOG_LEVEL >= 3 /* 3 because this file is usually very large */
    A.ParallelWriteMM(getmatfname("A.mtx").c_str(), true);
    #endif

    MPI_Barrier(comm);
    #endif
}

void ELBALogger::log_seed_matrix(CT<SharedSeeds>::PSpParMat& B)
{
    #if LOG_LEVEL >= 1
    size_t numseeds = B.getnnz();
    size_t numreads = B.getnrow();

    if (isroot) std::cout << "Overlap matrix B has " << numreads << " rows (readids), " << numreads << " columns (readids), and " << numseeds << " nonzeros (overlap seeds)\n" << std::endl;

    #if LOG_LEVEL >= 2
    B.ParallelWriteMM(getmatfname("B.mtx").c_str(), true, SharedSeeds::IOHandler());
    #endif

    MPI_Barrier(comm);
    #endif
}

void ELBALogger::log_overlap_matrix(CT<Overlap>::PSpParMat& R)
{
    #if LOG_LEVEL >= 2
    R.ParallelWriteMM(getmatfname("R.mtx").c_str(), true, Overlap::IOHandler());
    #endif
}

std::string ELBALogger::getmatfname(const std::string matname)
{
    std::ostringstream ss;
    ss << prefix << "." << matname;
    return ss.str();
}
