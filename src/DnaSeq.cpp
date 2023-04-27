#include "DnaSeq.hpp"
#include <cassert>
#include <cstring>
#include <vector>

DnaSeq::DnaSeq(char const *sequence, size_t len) : numbytes((len+3)/4), remain(4*numbytes - len), owns_memory(true)
{
    assert(remain < 4);

    memory = new uint8_t[numbytes];

    size_t b = 0;
    char const *sb = sequence;

    while (b < numbytes)
    {
        uint8_t byte = 0;
        int left = (b != numbytes-1? 4 : 4-remain);

        for (int i = 0; i < left; ++i)
        {
            uint8_t code = DnaSeq::getcharcode(sb[i]);
            uint8_t shift = code << (6 - (2*i));
            byte |= shift;
        }

        memory[b++] = byte;
        sb += 4;
    }
}

DnaSeq::DnaSeq(const DnaSeq& rhs) : numbytes(rhs.numbytes), remain(rhs.remain), owns_memory(true)
{
    memory = new uint8_t[numbytes];
    std::memcpy(memory, rhs.memory, numbytes);
}

std::string DnaSeq::ascii() const
{
    size_t len = size();
    uint8_t const *bb = memory;
    std::vector<char> s(len);

    for (size_t i = 0; i < len; ++i)
    {
        int code = (*bb >> (6 - (2*(i%4)))) & 3;
        s[i] = DnaSeq::getcodechar(code);

        if ((i+1) % 4 == 0)
            bb++;
    }

    return std::string(s.begin(), s.end());
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
