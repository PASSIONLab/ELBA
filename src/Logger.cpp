#include "Logger.hpp"
#include <iostream>
#include <numeric>

Logger::Logger(Grid commgrid) : logstream(new std::ostringstream()), commgrid(commgrid)
{
    comm = commgrid->GetWorld();
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &nprocs);
}

std::string Logger::readrangestr(size_t pos, size_t count)
{
    std::ostringstream ss;
    ss << "[" << pos << ".." << (pos+count-1) << "] (" << count << " reads)";
    return ss.str();
}

std::string Logger::rankstr(int proc)
{
    int nprocs = commgrid->GetSize();
    std::ostringstream ss;
    ss << "rank[" << proc+1 << "/" << nprocs << "]";
    return ss.str();
}

std::string Logger::rankstr()
{
    static bool initialized = false;
    static std::string name;

    if (!initialized)
    {
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();
        std::ostringstream ss;
        ss << "rank[" << myrank+1 << "/" << nprocs << "]";
        name = ss.str();
        initialized = true;
    }

    return name;
}

std::string Logger::prefix()
{
    static bool initialized = false;
    static std::string name;

    if (!initialized)
    {
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();
        std::ostringstream ss;
        ss << "rank[" << myrank+1 << "/" << nprocs << "] :: ";
        name = ss.str();
        initialized = true;
    }

    return name;
}

void Logger::Flush(std::ostringstream& ss)
{
    Flush(ss.str().c_str());
    ss.clear();
    ss.str("");
}

void Logger::Flush(std::ostringstream& ss, int rank)
{
    if (rank == myrank)
    {
        std::cout << ss.str() << std::endl;
    }

    ss.clear();
    ss.str("");

    MPI_Barrier(comm);
}

void Logger::Flush(char const *label)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    std::vector<int> recvcnt, displs;
    std::vector<char> recvbuf;

    std::string mylog = logstream->str();
    logstream.reset(new std::ostringstream());

    int sendcnt = mylog.size();
    int totrecv;

    if (!myrank) recvcnt.resize(nprocs);

    MPI_Gather(&sendcnt, 1, MPI_INT, recvcnt.data(), 1, MPI_INT, 0, comm);

    if (!myrank)
    {
        displs.resize(nprocs);
        displs.front() = 0;
        std::partial_sum(recvcnt.begin(), recvcnt.end()-1, displs.begin()+1);
        recvbuf.resize(recvcnt.back() + displs.back());
    }

    MPI_Gatherv(mylog.c_str(), sendcnt, MPI_CHAR, recvbuf.data(), recvcnt.data(), displs.data(), MPI_CHAR, 0, comm);

    if (!myrank)
    {
        std::string slabel(label);
        std::string banner;
        banner.assign(slabel.size(), '=');
        std::cout << slabel << "\n" << banner << "\n" << std::endl;

        char const *buf = recvbuf.data();

        for (int i = 0; i < nprocs; ++i)
        {
            std::cout << "rank[" << i+1 << "/" << nprocs << "] :: " << std::string(buf + displs[i], recvcnt[i]) << "\n";
        }
        std::cout << std::endl;
    }

    MPI_Barrier(comm);
}

void LogAll(const std::string mylog, Grid commgrid)
{
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    std::vector<int> recvcnt, displs;
    std::vector<char> recvbuf;
    int sendcnt = mylog.size();
    int totrecv;

    if (!myrank) recvcnt.resize(nprocs);

    MPI_Gather(&sendcnt, 1, MPI_INT, recvcnt.data(), 1, MPI_INT, 0, comm);

    if (!myrank)
    {
        displs.resize(nprocs);
        std::exclusive_scan(recvcnt.begin(), recvcnt.end(), displs.begin(), static_cast<int>(0));
        recvbuf.resize(recvcnt.back() + displs.back());
    }

    MPI_Gatherv(mylog.c_str(), sendcnt, MPI_CHAR, recvbuf.data(), recvcnt.data(), displs.data(), MPI_CHAR, 0, comm);

    if (!myrank)
    {
        char const *buf = recvbuf.data();

        for (int i = 0; i < nprocs; ++i)
        {
            std::string message(buf + displs[i], recvcnt[i]);
            std::cout << "rank[" << i+1 << "/" << nprocs << "] :: " << message << std::endl;
        }
        std::cout << std::endl;
    }

    MPI_Barrier(comm);
}

std::string ProcessorName(Grid commgrid)
{
    static bool initialized = false;
    static std::string name;

    if (!initialized)
    {
        int myrank = commgrid->GetRank();
        int nprocs = commgrid->GetSize();
        int rowrank = commgrid->GetRankInProcCol();
        int colrank = commgrid->GetRankInProcRow();
        std::ostringstream ss;
        ss  << "processor[" << myrank+1 << "/" << nprocs << "]...grid[" << rowrank+1 <<"," << colrank+1 << "])";
        name = ss.str();
        initialized = true;
    }

    return name;
}
