#include "DnaBuffer.hpp"

size_t DnaBuffer::computebufsize(const std::vector<size_t>& seqlens)
{
    auto bytecounter = [](size_t sum, size_t len) { return sum + DnaSeq::bytesneeded(len); };
    return std::accumulate(seqlens.cbegin(), seqlens.cend(), static_cast<size_t>(0), bytecounter);
}

void DnaBuffer::push_back(char const *s, size_t len) {}

// class DnaBuffer
// {
    // public:
    // DnaBuffer(size_t totbases, size_t totseqs);

    // size_t getnumseqs() const { return sequences.size(); }
    // size_t getbufsize() const { return buffer.size(); }

    // uint8_t* pushbufhead(size_t seqlen);

    // void push_back(char const *sequence, size_t len)
    // {
            // sequences.emplace_back(sequence, len, pushbufhead(len));
        // }

    // const DnaSeq& operator[](size_t i) { return sequences[i]; }

    // private:
    // size_t totbases; [> total number of nucleotides allowed <]
    // size_t totseqs; [> total number of sequences allowed <]
    // std::vector<uint8_t> buffer;
    // std::vector<DnaSeq> sequences;
// };
