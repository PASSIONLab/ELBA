#include "DnaBuffer.hpp"
#include "Logger.hpp"
#include <cassert>

size_t DnaBuffer::computebufsize(const std::vector<size_t>& seqlens)
{
    auto bytecounter = [](size_t sum, size_t len) { return sum + DnaSeq::bytesneeded(len); };
    return std::accumulate(seqlens.cbegin(), seqlens.cend(), static_cast<size_t>(0), bytecounter);
}

void DnaBuffer::push_back(char const *s, size_t len)
{
    size_t nbytes = DnaSeq::bytesneeded(len);
    assert(bufhead + nbytes <= bufsize);
    sequences.emplace_back(s, len, &buf[bufhead]);
    assert(nbytes == sequences.back().numbytes());
    bufhead += nbytes;
}

size_t DnaBuffer::getrangebufsize(size_t start, size_t count) const
{
    const auto& startseq = sequences[start];
    const auto& endseq = sequences[start+count-1];

    size_t startoffset = startseq.data() - &buf[0];
    size_t endoffset = (endseq.data() + endseq.numbytes()) - &buf[0];

    return startoffset-endoffset+1;
}
