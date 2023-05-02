#include "DnaSeq.hpp"
#include <cassert>
#include <cstring>
#include <vector>
#include <cmath>

void DnaSeq::compress(char const *s)
{
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
