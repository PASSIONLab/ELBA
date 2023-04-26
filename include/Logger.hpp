#ifndef LOGGER_H_
#define LOGGER_H_

#include "common.h"

void LogAll(const std::string mylog, Grid commgrid);
std::string ProcessorName(Grid commgrid);

class Logger
{
    std::unique_ptr<std::ostringstream> logstream, rootstream;
    Grid commgrid;
    int myrank, nprocs;
    MPI_Comm comm;

    std::string prefix();

public:
    Logger(Grid commgrid);
    void Flush(char const *label);
    void Flush(std::ostringstream& ss);
    void Flush(std::ostringstream& ss, int rank);
    std::ostringstream& operator()() { return *logstream; }
    std::string rankstr();
    std::string rankstr(int proc);
    static std::string readrangestr(size_t pos, size_t count);
};


#endif
