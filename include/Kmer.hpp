#ifndef KMER_H_
#define KMER_H_

#include <string>
#include <array>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdint>
#include <cstring>
#include "common.h"
#include "compiletime.h"
#include "HyperLogLog.hpp"
#include "DnaSeq.hpp"

template <int NLONGS>
class Kmer
{
public:

    static_assert(NLONGS != 0);

    static constexpr int NBYTES = 8 * NLONGS;

    typedef std::array<uint64_t, NLONGS> MERARR;
    typedef std::array<uint8_t,  NBYTES> BYTEARR;

    Kmer();
    Kmer(const DnaSeq& s);
    Kmer(char const *s);
    Kmer(const void *mem);
    Kmer(const Kmer& o);

    Kmer& operator=(const Kmer& o);

    std::string GetString() const;

    bool operator<(const Kmer& o) const;
    bool operator==(const Kmer& o) const;
    bool operator!=(const Kmer& o) const;

    Kmer GetExtension(int code) const;
    Kmer GetTwin() const;
    Kmer GetRep() const;

    uint64_t GetHash() const;
    const void* GetBytes() const { return reinterpret_cast<const void*>(longs.data()); }

    void CopyDataInto(void *mem) const { std::memcpy(mem, longs.data(), NBYTES); }
    void CopyDataFrom(const void *mem) { std::memcpy(longs.data(), mem, NBYTES); }

    static std::vector<Kmer> GetKmers(const DnaSeq& s);
    static std::vector<Kmer> GetRepKmers(const DnaSeq& s);

    template <int N>
    friend std::ostream& operator<<(std::ostream& os, const Kmer<N>& kmer);

private:

    union { MERARR  longs;
            BYTEARR bytes; };

    void set_kmer(const DnaSeq& s);
    void set_kmer(char const *s, bool const revcomp = false);
};

template <int NLONGS>
std::ostream& operator<<(std::ostream& os, const Kmer<NLONGS>& kmer)
{
    os << KMER_SIZE << "-mer(" << kmer.GetString() << ")";
    return os;
}

namespace std
{
    template <int NLONGS> struct hash<Kmer<NLONGS>>
    {
        size_t operator()(const Kmer<NLONGS>& kmer) const
        {
            auto myhash = kmer.GetHash();
            return myhash;
        }
    };

    template <int NLONGS> struct less<Kmer<NLONGS>>
    {
        bool operator()(const Kmer<NLONGS>& k1, const Kmer<NLONGS>& k2) const
        {
            return k1 < k2;
        }
    };
}

#include "Kmer.cpp"

using TKmer = typename std::conditional<(KMER_SIZE <= 32), Kmer<1>,
              typename std::conditional<(KMER_SIZE <= 64), Kmer<2>,
              typename std::conditional<(KMER_SIZE <= 96), Kmer<3>, Kmer<0>>::type>::type>::type;

#endif
