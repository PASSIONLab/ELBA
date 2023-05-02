#ifndef DNABUFFER_H_
#define DNABUFFER_H_

#include "DnaSeq.hpp"
#include <memory>
#include <numeric>

class DnaBuffer
{
public:
    DnaBuffer(size_t bufsize) : bufhead(0), bufsize(bufsize), buf(new uint8_t[bufsize]) {}
    DnaBuffer(size_t bufsize, size_t numreads, uint8_t *buf, const size_t *readlens); /* TODO */

    void push_back(char const *s, size_t len);
    size_t size() const { return sequences.size(); }
    size_t getbufsize() const { return bufsize; }
    size_t getrangebufsize(size_t start, size_t count) const;
    const uint8_t* getbufoffset(size_t i) const { return sequences[i].data(); }
    const DnaSeq& operator[](size_t i) const { return sequences[i]; }

    static size_t computebufsize(const std::vector<size_t>& seqlens);

    ~DnaBuffer() { delete[] buf; }

private:
    size_t bufhead;
    const size_t bufsize;
    uint8_t *buf;
    std::vector<DnaSeq> sequences;
};

#endif
