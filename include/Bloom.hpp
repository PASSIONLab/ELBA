#ifndef BLOOM_H_
#define BLOOM_H_

#include <cstdint>
#include <cstddef>

class Bloom
{
public:
    Bloom(int64_t entries, double error);
    ~Bloom();

    bool Check(const void *buffer, size_t len);
    bool Add(const void *buffer, size_t len);

    int64_t entries;
    int64_t bits;
    int64_t bytes;
    double error;
    int hashes;

private:
    double bpe;
    unsigned char *bf;
    bool ready;

    bool bloom_check_add(const void *buffer, size_t len, bool add);
};

#endif
