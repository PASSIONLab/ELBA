#include "DnaBuffer.hpp"
#include "Logger.hpp"
#include <cassert>

// void DnaBufferContainer::assign(const std::vector<DnaSeq>& sequences)
// {
    // assert(sequences.size() == numseqs);

    // size_t bufhead = 0;

    // for (size_t i = 0; i < numseqs; ++i)
    // {
        // sequences[i].copyto(readlens + i, buf + bufhead);
        // bufhead += sequences[i].numbytes();
    // }
// }

// DnaBuffer::DnaBuffer(DnaBufferContainer& container) : bufhead(0), bufsize(container.bufsize), buf(container.buf)
// {
    // sequences.reserve(container.numseqs);

    // for (size_t i = 0; i < container.numseqs; ++i)
    // {
        // sequences.emplace_back(container.readlens[i], buf + bufhead);
        // bufhead += sequences.back().numbytes();
    // }

    // container.releasebuffer();
// }

size_t DnaBuffer::computebufsize(const std::vector<size_t>& seqlens)
{
    auto bytecounter = [](size_t sum, size_t len) { return sum + DnaSeq::bytesneeded(len); };
    return std::accumulate(seqlens.cbegin(), seqlens.cend(), static_cast<size_t>(0), bytecounter);
}

void DnaBuffer::push_back(char const *s, size_t len)
{
    size_t nbytes = DnaSeq::bytesneeded(len);
    assert(bufhead + nbytes <= bufsize);
    sequences.emplace_back(s, len, buf + bufhead);
    assert(nbytes == sequences.back().numbytes());
    bufhead += nbytes;
}

size_t DnaBuffer::getrangebufsize(size_t start, size_t count) const
{
    if (start + count == 0) return 0;
    size_t end = start+count-1;
    const uint8_t* startmem = sequences[start].data();
    const uint8_t* endmem = sequences[end].data() + sequences[end].numbytes();
    return (endmem-startmem);
}
