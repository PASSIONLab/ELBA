#include "DnaSeq.hpp"
#include <cassert>
#include <cstring>
#include <vector>
#include <cmath>

/*
 * For a sequence of length l, floor((l+3)/4) bytes are needed to encode it. Let
 * L = l[0] + l[1] + ... + l[N-1], where l[i] is the length of the ith sequence, N is
 * the total number of sequences, and hence L is the total sum of all the sequence
 * lengths. Then the total number of bytes needed for the write buffer is
 *
 * Sum{0 <= i <= N-1}[floor((l[i]+3)/4)] <= (1/4) * Sum{0 <= i <= N-1}[l[i] + 4]
 *                                       <= (1/4) * (L + 4N)
 *                                        = (L/4) + N
 */
DnaBuffer::DnaBuffer(size_t totbases, size_t totseqs) : totbases(totbases), totseqs(totseqs), numseqs(0)
{
    size_t bufmem = std::ceil((totbases/4.0)+totseqs);
    buffer.reserve(bufmem);
}

uint8_t* DnaBuffer::pushbufhead(size_t seqlen)
{
    assert(numseqs++ < totseqs);
    size_t bufsize = getbufsize();
    uint8_t *head = buffer.data() + bufsize;
    buffer.resize(bufsize + DnaSeq::bytesneeded(seqlen));
    return head;
}

DnaSeq::DnaSeq(char const *s, size_t len, DnaBuffer& extbuf) : len(len), ownsmem(false)
{
    memory = extbuf.pushbufhead(len);

    const size_t nbytes = numbytes();
    const int remain = remainder();
    char const *p = s;
    size_t b = 0;

    while (b < nbytes)
    {
        uint8_t byte = 0;
        int left = (b != nbytes-1? 4 : 4-remain);

        for (int i = 0; i < left; ++i)
        {
            uint8_t code = DnaSeq::getcharcode(p[i]);
            uint8_t shift = code << (6 - (2*i));
            byte |= shift;
        }

        memory[b++] = byte;
        p += 4;
    }
}

DnaSeq::DnaSeq(char const *s, size_t len) : len(len), ownsmem(true)
{
    memory = new uint8_t[numbytes()];

    const size_t nbytes = numbytes();
    const int remain = remainder();
    char const *p = s;
    size_t b = 0;

    while (b < nbytes)
    {
        uint8_t byte = 0;
        int left = (b != nbytes-1? 4 : 4-remain);

        for (int i = 0; i < left; ++i)
        {
            uint8_t code = DnaSeq::getcharcode(p[i]);
            uint8_t shift = code << (6 - (2*i));
            byte |= shift;
        }

        memory[b++] = byte;
        p += 4;
    }
}

DnaSeq::DnaSeq(const DnaSeq& rhs) : len(rhs.len), ownsmem(true)
{
    memory = new uint8_t[numbytes()];
    std::memcpy(memory, rhs.memory, numbytes());
}

std::string DnaSeq::ascii() const
{
    size_t len = size();
    uint8_t const *p = memory;
    std::string s(len, '\0');

    for (size_t i = 0; i < len; ++i)
    {
        int code = (*p >> (6 - (2*(i%4)))) & 3;
        s[i] = DnaSeq::getcodechar(code);

        if ((i+1) % 4 == 0)
            p++;
    }
    return s;
}


int DnaSeq::operator[](size_t i) const
{
    uint8_t byte = memory[i/4];
    int shift = 6 - (2 * (i%4));
    int code = (byte >> shift)&3;
    return code;
}

bool DnaSeq::operator==(const DnaSeq& rhs)
{
    if (size() != rhs.size())
        return false;

    for (size_t i = 0; i < rhs.size(); ++i)
        if ((*this)[i] != rhs[i])
            return false;

    return true;
}

bool DnaSeq::operator<(const DnaSeq& rhs)
{
    size_t len = std::min(size(), rhs.size());
    size_t i;

    for (i = 0; (*this)[i] == rhs[i] && i < len; ++i);

    if (i == len) return false;

    return ((*this)[i] < rhs[i]);
}
