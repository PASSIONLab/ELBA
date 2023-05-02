#ifndef READ_OVERLAP_H_
#define READ_OVERLAP_H_

#include "KmerOps.hpp"
#include <iostream>

struct ReadOverlap
{
    int sfx, sfxT, dir, dirT;
    int score;
    PosInRead begQs[2], begTs[2];
    PosInRead b[2], e[2], l[2], coords[2];
    int sfxpath[4];
    int count;
    bool transpose, passed, rc;

    void SetPathInf();

    ReadOverlap();
    ReadOverlap(int count);
    ReadOverlap(const ReadOverlap& rhs);

    bool is_invalid() const;

    bool arrows(int& t, int& h) const;

    ReadOverlap operator+(const ReadOverlap& b)
    {
        ReadOverlap o = b;
        return o;
    }

    int intplus(int a, int b);

    friend std::ostream& operator<<(std::ostream& os, const ReadOverlap& o)
    {
        os << o.begQs[0] << "\t" << o.begTs[0] << "\t" << o.begQs[1] << "\t" << o.begTs[1];
        return os;
    }
};

struct OverlapHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& o, uint64_t row, uint64_t col)
    {
        os << o;
    }
};

#endif
