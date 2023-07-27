#include "HyperLogLog.hpp"
#include "HashFuncs.hpp"
#include <cmath>
#include <cassert>

/* reference: Aydin Buluc modified version of HyperLogLog written by Hideaki Ohno */

#define HASHBITS (64)
#define SIGNBIT (1ULL << (HASHBITS-1))
#define BIG32 ((double)(1ULL << 32))

static inline uint8_t rho(uint64_t hash, uint8_t bits)
{
    uint8_t v = 1;

    while (v <= bits && !(hash & SIGNBIT))
    {
        v++;
        hash <<= 1;
    }

    return v;
}

HyperLogLog::HyperLogLog(uint8_t bits) : bits(bits), size(1ULL << bits), registers(size+1, 0)
{
    double alpha;

    switch (size)
    {
        case 16: alpha = 0.673; break;
        case 32: alpha = 0.697; break;
        case 64: alpha = 0.709; break;
        default: alpha = 0.7213 / (1.0 + (1.079/size));
    }

    alpha_mm = alpha * size * size;
}

void HyperLogLog::add(char const *s, size_t len)
{
    uint64_t hashval;
    murmurhash3_64(s, len, &hashval);

    uint32_t index = hashval >> (HASHBITS - bits);
    uint8_t rank = rho((hashval << bits), HASHBITS - bits);

    if (rank > registers[index])
        registers[index] = rank;
}

double HyperLogLog::estimate() const
{
    double est;
    double sum = 0.0;

    for (uint32_t i = 0; i < size; ++i)
        sum += 1.0 / (1 << registers[i]);

    est = alpha_mm / sum;

    if (est <= 2.5 * size)
    {
        uint32_t zeros = 0;

        for (uint32_t i = 0; i < size; ++i)
            zeros += (registers[i] == 0);

        if (zeros)
        {
            est = size * log((double)size / zeros);
        }
    }

    return est;
}

void HyperLogLog::merge(const HyperLogLog& rhs)
{
    assert(bits == rhs.bits);

    for (uint32_t i = 0; i < size; ++i)
    {
        registers[i] = std::max(registers[i], rhs.registers[i]);
    }
}

HyperLogLog& HyperLogLog::parallelmerge(MPI_Comm comm)
{
    MPI_ALLREDUCE(MPI_IN_PLACE, registers.data(), static_cast<MPI_Count_type>(size+1), MPI_UINT8_T, MPI_MAX, comm);
    return *this;
}
