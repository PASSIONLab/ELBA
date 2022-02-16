#ifndef READ_OVERLAPS_H_
#define READ_OVERLAPS_H_

#include <limits>
#include <cassert>
#include <iostream>
#include "kmer/CommonKmers.hpp"

static constexpr int MAX_INT = std::numeric_limits<int>::max();
static constexpr int XBOUND = 150;

struct ReadOverlap
{
    int rc;
    int b[2], e[2], l[2];
    int sfx[4];
    int valid;
    int overlap;

    ReadOverlap() : valid(1), rc(0), overlap(0), b{}, e{}, l{}
    {
        for (int i = 0; i < 4; ++i) sfx[i] = MAX_INT;
    }

    ReadOverlap(const ReadOverlap& rhs) : rc(rhs.rc), valid(rhs.valid), overlap(rhs.overlap)
    {
        int i;

        for (i = 0; i < 4; ++i)
            sfx[i] = rhs.sfx[i];

        for (i = 0; i < 2; ++i) {
            b[i] = rhs.b[i];
            e[i] = rhs.e[i];
            l[i] = rhs.l[i];
        }
    }

    ReadOverlap(const dibella::CommonKmers& cks) : valid(1), overlap(0)
    {
        b[0] = cks.first.first;  b[1] = cks.second.first;
        e[0] = cks.first.second; e[1] = cks.second.second;
        l[0] = cks.lenv;         l[1] = cks.lenh;

        rc = cks.rc;
        overlap = cks.overlap;

        refix();
    }

    void refix(int sdove = 0)
    {
        if (!valid) return;

        for (int i = 0; i < 4; ++i) sfx[i] = MAX_INT;

        if (l[0] - e[0] < XBOUND && b[1] < XBOUND)
            sfx[(rc) ? (0) : (1+sdove)] = l[1] - e[1];
        else if (l[1] - e[1] < XBOUND && b[0] < XBOUND)
            sfx[(rc) ? (3) : (2-sdove)] = b[1];
        else
            valid = 0;
    }

    int direction() const
    {
        for (int i = 0; i < 4; ++i) { if (sfx[i] < MAX_INT) return i; }
        return -1;
    }

    int getsuffix(int& dir) const
    {
        dir = direction();

        return (dir != -1) ? sfx[dir] : MAX_INT;
    }

    friend bool operator==(const ReadOverlap& lhs, const ReadOverlap& rhs)
    {
        for (int i = 0; i < 4; ++i)
            if (lhs.sfx[i] != rhs.sfx[i]) return false;
        return true;
    }

    ReadOverlap operator+(const ReadOverlap& b)
    {
        ReadOverlap myobj;
        myobj = b;
        return myobj;
    }
};

struct ReadOverlapHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, uint64_t row, uint64_t col)
    {
        os << e.direction() << "\t" << e.sfx[0] << "\t" << e.sfx[1] << "\t" << e.sfx[2] << "\t" << e.sfx[3];
    }
};

struct ReadOverlapMMHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, uint64_t row, uint64_t col)
    {
        int dir, suf;
        suf = e.getsuffix(dir);
        os << dir << "\t" << suf;
    }
};

#endif
