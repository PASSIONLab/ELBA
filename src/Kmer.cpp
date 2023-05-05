#include "Kmer.hpp"
#include "HashFuncs.hpp"
#include "DnaSeq.hpp"
#include <cstring>
#include <cassert>
#include <algorithm>
#include <limits>

static uint64_t tetramer_twin(const uint8_t code)
{
    static const uint8_t tetramer_lookup_code[256] =
    {
        0xff,0xbf,0x7f,0x3f,0xef,0xaf,0x6f,0x2f,0xdf,0x9f,0x5f,0x1f,0xcf,0x8f,0x4f,0xf,
        0xfb,0xbb,0x7b,0x3b,0xeb,0xab,0x6b,0x2b,0xdb,0x9b,0x5b,0x1b,0xcb,0x8b,0x4b,0xb,
        0xf7,0xb7,0x77,0x37,0xe7,0xa7,0x67,0x27,0xd7,0x97,0x57,0x17,0xc7,0x87,0x47,0x7,
        0xf3,0xb3,0x73,0x33,0xe3,0xa3,0x63,0x23,0xd3,0x93,0x53,0x13,0xc3,0x83,0x43,0x3,
        0xfe,0xbe,0x7e,0x3e,0xee,0xae,0x6e,0x2e,0xde,0x9e,0x5e,0x1e,0xce,0x8e,0x4e,0xe,
        0xfa,0xba,0x7a,0x3a,0xea,0xaa,0x6a,0x2a,0xda,0x9a,0x5a,0x1a,0xca,0x8a,0x4a,0xa,
        0xf6,0xb6,0x76,0x36,0xe6,0xa6,0x66,0x26,0xd6,0x96,0x56,0x16,0xc6,0x86,0x46,0x6,
        0xf2,0xb2,0x72,0x32,0xe2,0xa2,0x62,0x22,0xd2,0x92,0x52,0x12,0xc2,0x82,0x42,0x2,
        0xfd,0xbd,0x7d,0x3d,0xed,0xad,0x6d,0x2d,0xdd,0x9d,0x5d,0x1d,0xcd,0x8d,0x4d,0xd,
        0xf9,0xb9,0x79,0x39,0xe9,0xa9,0x69,0x29,0xd9,0x99,0x59,0x19,0xc9,0x89,0x49,0x9,
        0xf5,0xb5,0x75,0x35,0xe5,0xa5,0x65,0x25,0xd5,0x95,0x55,0x15,0xc5,0x85,0x45,0x5,
        0xf1,0xb1,0x71,0x31,0xe1,0xa1,0x61,0x21,0xd1,0x91,0x51,0x11,0xc1,0x81,0x41,0x1,
        0xfc,0xbc,0x7c,0x3c,0xec,0xac,0x6c,0x2c,0xdc,0x9c,0x5c,0x1c,0xcc,0x8c,0x4c,0xc,
        0xf8,0xb8,0x78,0x38,0xe8,0xa8,0x68,0x28,0xd8,0x98,0x58,0x18,0xc8,0x88,0x48,0x8,
        0xf4,0xb4,0x74,0x34,0xe4,0xa4,0x64,0x24,0xd4,0x94,0x54,0x14,0xc4,0x84,0x44,0x4,
        0xf0,0xb0,0x70,0x30,0xe0,0xa0,0x60,0x20,0xd0,0x90,0x50,0x10,0xc0,0x80,0x40,0x0
    };

    return static_cast<uint64_t>(tetramer_lookup_code[code]);
}

template <int NLONGS>
Kmer<NLONGS>::Kmer() : longs{} {}

template <int NLONGS>
Kmer<NLONGS>::Kmer(const DnaSeq& s) : Kmer() { set_kmer(s); }

template <int NLONGS>
Kmer<NLONGS>::Kmer(char const *s) : Kmer() { set_kmer(s); }

template <int NLONGS>
Kmer<NLONGS>::Kmer(const Kmer& o) : longs(o.longs) {}

template <int NLONGS>
Kmer<NLONGS>::Kmer(const void *mem) : Kmer() { CopyDataFrom(mem); }

template <int NLONGS>
std::string Kmer<NLONGS>::GetString() const
{
    std::string s(KMER_SIZE, '\0');

    int i, j, l;

    for (i = 0; i < KMER_SIZE; ++i)
    {
        j = i % 32;
        l = i / 32;

        s[i] = "ACGT"[(longs[l] >> (2 * (31 - j)))&3];
    }

    return s;
}

template <int NLONGS>
void Kmer<NLONGS>::set_kmer(const DnaSeq& s)
{
    int i, j, l, idx;
    uint64_t code;

    /*
     * Warning: set_kmer assumes that longs/bytes have already
     * been completely zeroed out.
     */

    for (i = 0; i < KMER_SIZE; ++i)
    {
        j = i % 32;
        l = i / 32;

        code = static_cast<uint64_t>(s[i]);

        longs[l] |= (code << (2 * (31 - j)));
    }
}

template <int NLONGS>
void Kmer<NLONGS>::set_kmer(char const *s, bool const revcomp)
{
    int i, j, l, idx;
    uint64_t code;

    /*
     * Warning: set_kmer assumes that longs/bytes have already
     * been completely zeroed out.
     */

    for (i = 0; i < KMER_SIZE; ++i)
    {
        j = i % 32;
        l = i / 32;

        idx = revcomp? KMER_SIZE - i - 1 : i;
        code = static_cast<uint64_t>(DnaSeq::getcharcode(s[idx]));

        longs[l] |= ((revcomp? 3 - code : code) << (2 * (31 - j)));
    }
}
template <int NLONGS>
Kmer<NLONGS>& Kmer<NLONGS>::operator=(Kmer o)
{
    std::swap(longs, o.longs);
    return *this;
}

template <int NLONGS>
bool Kmer<NLONGS>::operator<(const Kmer& o) const
{
    for (int i = 0; i < NLONGS; ++i)
    {
        if (longs[i] < o.longs[i])
            return true;

        if (longs[i] > o.longs[i])
            return false;
    }

    return false;
}

template <int NLONGS>
bool Kmer<NLONGS>::operator==(const Kmer& o) const
{
    for (int i = 0; i < NLONGS; ++i)
        if (longs[i] != o.longs[i])
            return false;

    return true;
}

template <int NLONGS>
bool Kmer<NLONGS>::operator!=(const Kmer& o) const
{
    return !(*this == o);
}

template <int NLONGS>
Kmer<NLONGS> Kmer<NLONGS>::GetExtension(int code) const
{
    Kmer ext;

    ext.longs[0] = longs[0] << 2;

    for (uint64_t i = 1; i < NLONGS; ++i)
    {
        ext.longs[i-1] |= ((longs[i] >> 62) & 0x3);
        ext.longs[i] = longs[i] << 2;
    }

    ext.longs[NLONGS-1] |= (static_cast<uint64_t>(code) << (2 * (32 - (KMER_SIZE%32))));

    return ext;
}

template <int NLONGS>
Kmer<NLONGS> Kmer<NLONGS>::GetTwin() const
{
    Kmer twin;

    /* unroll */
    for (int l = 0; l < NLONGS; ++l)
    {
        uint64_t longmer = longs[l];

        /* unroll */
        for (uint64_t i = 0; i < 64; i += 8)
        {
            uint8_t bytemer = (longmer >> i) & 0xff;
            uint64_t revcomp_bytemer = tetramer_twin(bytemer);
            twin.longs[NLONGS-1-l] |= (revcomp_bytemer << (56 - i));
        }
    }

    uint64_t shift = KMER_SIZE % 32? 2 * (32 - (KMER_SIZE % 32)) : 0ULL;
    uint64_t mask = KMER_SIZE % 32? ((1ULL << shift) - 1) << (64 - shift) : 0ULL;

    twin.longs[0] <<= shift;

    for (uint64_t i = 1; i < NLONGS; ++i)
    {
        twin.longs[i-1] |= (twin.longs[i] & mask) >> (64 - shift);
        twin.longs[i] <<= shift;
    }

    return twin;
}

template <int NLONGS>
Kmer<NLONGS> Kmer<NLONGS>::GetRep() const
{
    Kmer twin = GetTwin();
    return twin < *this? twin : *this;
}

template <int NLONGS>
uint64_t Kmer<NLONGS>::GetHash() const
{
    uint64_t h;
    murmurhash3_64(longs.data(), NBYTES, &h);
    return h;
}

template <int NLONGS>
std::vector<Kmer<NLONGS>> Kmer<NLONGS>::GetKmers(const DnaSeq& s)
{
    int l = s.size();
    int num_kmers = l - KMER_SIZE + 1;

    if (num_kmers <= 0) return std::vector<Kmer>();

    std::vector<Kmer> kmers;

    kmers.reserve(num_kmers);
    kmers.emplace_back(s);

    for (int i = 1; i < num_kmers; ++i)
    {
        kmers.push_back(kmers.back().GetExtension(s[i+KMER_SIZE-1]));
    }

    return kmers;
}

template <int NLONGS>
std::vector<Kmer<NLONGS>> Kmer<NLONGS>::GetRepKmers(const DnaSeq& s)
{
    auto kmers = GetKmers(s);
    std::transform(kmers.begin(), kmers.end(), kmers.begin(), [](const Kmer& kmer) { return kmer.GetRep(); });
    return kmers;
}
