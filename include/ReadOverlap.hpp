#ifndef READ_OVERLAPS_H_
#define READ_OVERLAPS_H_

#include <limits>
#include <cassert>
#include <iostream>
#include "kmer/CommonKmers.hpp"

static constexpr int MAX_INT = std::numeric_limits<int>::max();
static constexpr int XBOUND = 150;

struct ReadOverlap;

struct OverlapPath
{
    int sfx[4];
    OverlapPath();
    OverlapPath(const ReadOverlap& e);
    void setinf();
};

struct ReadOverlap
{
    int sfx, dir, rc, overlap;
    int b[2], e[2], l[2];
    int transpose;

    int coords[2];

    ReadOverlap() : sfx(0), dir(-1), transpose(0) {}
    ReadOverlap(int sfx, int dir) : sfx(sfx), dir(dir), transpose(0)
    {
        if (sfx <= 0) dir = -1;
    }

    ReadOverlap(const OverlapPath& e) : transpose(0)
    {
        int c = 0;
        for (int i = 0; i < 4; ++i) {
            if (e.sfx[i] < MAX_INT) {
                sfx = e.sfx[i];
                dir = i;
                c++;
            }
        }

        if (c != 1) {
            sfx = 0;
            dir = -1;
        }
    }

    ReadOverlap(const ReadOverlap& lhs) : sfx(lhs.sfx), dir(lhs.dir), rc(lhs.rc), overlap(lhs.overlap), transpose(lhs.transpose)
    {
        for (int i = 0; i < 2; ++i) {
            b[i] = lhs.b[i];
            e[i] = lhs.e[i];
            l[i] = lhs.l[i];
            coords[i] = lhs.coords[i];
        }
    }

    ReadOverlap(const dibella::CommonKmers& cks)
    {
        b[0] = cks.first.first;  b[1] = cks.second.first;
        e[0] = cks.first.second; e[1] = cks.second.second;
        l[0] = cks.lenv;         l[1] = cks.lenh;

        rc = cks.rc;

        refix();
    }

    void refix(int sdove = 0)
    {
        if (dir == -1) return;

        if (sdove==1) transpose = 1;

        if (l[0] - e[0] < XBOUND && b[1] < XBOUND) {
            dir = (rc) ? (0) : (1+sdove);
            sfx = l[1] - e[1];
            overlap = e[0] - b[0];
        } else if (l[1] - e[1] < XBOUND && b[0] < XBOUND) {
            dir = (rc) ? (3) : (2-sdove);
            sfx = b[1];
            overlap = e[0] - b[0];
        } else {
            sfx = 0;
            dir = -1;
        }
    }

    bool isvalid() const { return (sfx > 0 && dir >= 0); }

    bool arrows(int& t, int& h) const
    {
        if (!isvalid())
            return false;

        t = (dir>>1)&1;
        h = dir&1;

        return true;
    }

    friend bool operator==(const ReadOverlap& lhs, const ReadOverlap& rhs)
    {
        return (lhs.sfx == rhs.sfx && lhs.dir == rhs.dir);
    }

    operator bool() const { return isvalid(); }

    ReadOverlap operator+(const ReadOverlap& b)
    {
        ReadOverlap myobj = b;
        return myobj;
    }
};

OverlapPath::OverlapPath() { setinf(); }
OverlapPath::OverlapPath(const ReadOverlap& e)
{
    setinf();

    if (e.isvalid())
        sfx[e.dir] = e.sfx;
}

void OverlapPath::setinf() { for (int i = 0; i < 4; ++i) sfx[i] = MAX_INT; }

struct ReadOverlapExtraHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << e.dir << "\t" << e.transpose << "\t" << e.b[0] << "\t" << e.e[0] << "\t" << e.l[0] << "\t" << e.b[1] << "\t" << e.e[1] << "\t" << e.l[1];
    }
    
};

struct ReadOverlapHandler
{
    ReadOverlap getNoNum(int64_t row, int64_t col) { return ReadOverlap(); }

    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << e.dir << "\t" << e.sfx;
    }

    template <typename c, typename t>
    void read(std::basic_istream<c,t>& is, const ReadOverlap& e, int64_t row, int64_t col)
    {
        int sfx, dir;
        is >> dir >> sfx;

        return ReadOverlap(sfx, dir);
    }

    void binaryfill(FILE *rFile, int64_t row, int64_t col, ReadOverlap& e) {}

    size_t entrylength() { return sizeof(ReadOverlap); }
};

int intplus(int a, int b)
{
    return (a == MAX_INT || b == MAX_INT) ? MAX_INT : a + b;
}

struct Tupleize : unary_function<ReadOverlap, ReadOverlap>
{
    ReadOverlap operator() (ReadOverlap& e)
    {
        ReadOverlap out = e;
        switch (e.dir) {
            case 0:
                out.coords[0] = e.b[0] + 15;
                out.coords[1] = e.l[1] - e.b[1];
                break;
            case 3:
                out.coords[0] = e.e[0] - 15;
                out.coords[1] = e.l[1] - e.e[1];
                break;
            case 1:
                out.coords[0] = (e.transpose)? (e.l[0] - e.e[0] + 15) : (e.b[0] + 15);
                out.coords[1] = (e.transpose)? (e.l[1] - e.e[1]) : (e.b[1]);
                break;
            case 2:
                out.coords[0] = (e.transpose)? (e.l[0] - e.b[0] - 15) : (e.e[0] - 15);
                out.coords[1] = (e.transpose)? (e.l[1] - e.b[1]) : (e.e[1]);
                break;
            default:
                break;
        }
        return out;
    }
};
struct TupleHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << e.coords[0] << "\t" << e.coords[1];
    }
};


#endif
