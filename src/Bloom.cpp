#include "Bloom.hpp"
#include "HashFuncs.hpp"
#include <cassert>
#include <cmath>

Bloom::Bloom(int64_t entries, double error) : entries(entries), error(error), ready(false)
{
    assert(entries >= 1 && error > 0 && error < 1);

    double num = std::log(error);
    double denom = 0.480453013918201; // ln(2)^2

    bpe = -(num / denom);

    double dentries = (double)entries;
    bits = (int64_t)(dentries * bpe);

    bytes = (bits / 8) + !!(bits % 8);

    hashes = (int)std::ceil(0.693147180559945 * bpe);  // ln(2)

    bf = new unsigned char[bytes];

    std::fill_n(bf, bytes, static_cast<unsigned char>(0));

    ready = true;
}

Bloom::~Bloom()
{
    delete [] bf;
}

bool Bloom::Check(const void *buffer, size_t len)
{
    return bloom_check_add(buffer, len, false);
}

bool Bloom::Add(const void *buffer, size_t len)
{
    return bloom_check_add(buffer, len, true);
}

bool Bloom::bloom_check_add(const void *buffer, size_t len, bool add)
{
    assert(ready);

    int hits = 0;
    uint32_t a1 = murmurhash3(buffer, len, 0x9747b28c);
    uint32_t a2 = murmurhash3(buffer, len, a1);
    uint32_t b1 = murmurhash3(buffer, len, a2);
    uint32_t b2 = murmurhash3(buffer, len, b1);
    uint64_t a = (((uint64_t)a1)<<32) | ((uint64_t)a2);
    uint64_t b = (((uint64_t)b1)<<32) | ((uint64_t)b2);
    uint64_t x;
    uint64_t byte;
    uint32_t mask;
    uint32_t i;
    unsigned char c;

    for (i = 0; i < hashes; ++i)
    {
        x = (a + i*b) % bits;
        byte = x >> 3;
        c = bf[byte];
        mask = 1 << (x % 8);

        if (c & mask) hits++;
        else if (add) bf[byte] = c | mask;
    }

    return (hits == hashes);
}
