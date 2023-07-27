#ifndef HYPERLOGLOG_H
#define HYPERLOGLOG_H

#include "common.h"
#include <cstdint>
#include <mpi.h>

class HyperLogLog
{
private:
    uint8_t bits;    /// register bit width
    uint32_t size;   /// register size
    double alpha_mm; /// alpha * m^2
    std::vector<uint8_t> registers;

public:
    HyperLogLog(uint8_t bits = 12);
    void add(const char* s, size_t len);
    void add(const std::string& s) { add(s.c_str(), s.size()); }
    double estimate() const;
    void merge(const HyperLogLog& rhs);
    HyperLogLog& parallelmerge(MPI_Comm comm);
};

#endif
