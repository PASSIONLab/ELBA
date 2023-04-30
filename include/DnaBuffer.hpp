#ifndef DNABUFFER_H_
#define DNABUFFER_H_

#include "DnaSeq.hpp"
#include <memory>
#include <numeric>

class DnaBuffer
{
public:
    DnaBuffer(size_t bufsize) : bufsize(bufsize), bufhead(0), buf(new uint8_t[bufsize]) {}

    void push_back(char const *s, size_t len);
    size_t size() const { return sequences.size(); }
    size_t getbufsize() const { return bufsize; }
    const DnaSeq& operator[](size_t i) { return sequences[i]; }
    static size_t computebufsize(const std::vector<size_t>& seqlens);

private:
    size_t bufsize, bufhead;
    std::unique_ptr<uint8_t[]> buf;
    std::vector<DnaSeq> sequences;
};

#endif
