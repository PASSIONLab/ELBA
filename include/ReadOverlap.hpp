#ifndef READ_OVERLAPS_H_
#define READ_OVERLAPS_H_

#include <limits>
#include <cassert>
#include <iostream>
#include "kmer/CommonKmers.hpp"

using namespace dibella;

static constexpr int MAX_INT = std::numeric_limits<int>::max();
extern int xdrop;

struct ReadOverlap;

struct OverlapPath
{
    int sfx[4];
    bool transpose;
    OverlapPath() : transpose(false) { setinf(); }
    OverlapPath(const ReadOverlap& e);
    void setinf();
};

struct ReadOverlap
{
    int sfx, sfxT, dir, dirT;
    int b[2], e[2], l[2], coords[2];
    bool transpose, rc;

    ReadOverlap() : sfx(0), dir(-1), transpose(false) {}
    ReadOverlap(const ReadOverlap& rhs)
        : sfx(rhs.sfx), sfxT(rhs.sfxT), dir(rhs.dir), dirT(rhs.dirT), transpose(rhs.transpose), rc(rhs.rc)
    {
        b[0] = rhs.b[0]; b[1] = rhs.b[1];
        e[0] = rhs.e[0]; e[1] = rhs.e[1];
        l[0] = rhs.l[0]; l[1] = rhs.l[1];

        coords[0] = rhs.coords[0];
        coords[1] = rhs.coords[1];
    }

    ReadOverlap(const CommonKmers& cks) : transpose(false)
    {
        b[0] = cks.first.first;  b[1] = cks.second.first;
        e[0] = cks.first.second; e[1] = cks.second.second;
        l[0] = cks.lenv;         l[1] = cks.lenh;

        rc = cks.rc;

        sfx  = cks.sfx;
        sfxT = cks.sfxT;
        dir  = cks.dir;
        dirT = cks.dirT;
    }

    ReadOverlap(const OverlapPath& e) : transpose(e.transpose)
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

    bool is_invalid() const { return (dir == -1); }

    bool arrows(int& t, int& h) const
    {
        if (is_invalid())
            return false;

        t = (dir>>1)&1;
        h = dir&1;

        return true;
    }

    operator bool() const { return true; } /* for T = R in transitive reduction */

    //operator int64_t() const { return (is_invalid())? 0 : 1; }

    friend bool operator==(const ReadOverlap& lhs, const ReadOverlap& rhs)
    {
        return (lhs.sfx == rhs.sfx && lhs.dir == rhs.dir);
    }

    ReadOverlap operator+(const ReadOverlap& b)
    {
        ReadOverlap myobj = b;
        return myobj;
    }
};

OverlapPath::OverlapPath(const ReadOverlap& e) : transpose(e.transpose)
{
    setinf();

    if (!e.is_invalid())
        sfx[e.dir] = e.sfx;
}

void OverlapPath::setinf() { sfx[0] = sfx[1] = sfx[2] = sfx[3] = MAX_INT; }

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
                out.coords[0] = e.b[0] + xdrop;
                out.coords[1] = e.l[1] - e.b[1];
                break;
            case 1:
                out.coords[0] = (e.transpose)? (e.l[0] - e.e[0] + xdrop) : (e.b[0] + xdrop);
                out.coords[1] = (e.transpose)? (e.l[1] - e.e[1]) : (e.b[1]);
                break;
            case 2:
                out.coords[0] = (e.transpose)? (e.l[0] - e.b[0] - xdrop) : (e.e[0] - xdrop);
                out.coords[1] = (e.transpose)? (e.l[1] - e.b[1]) : (e.e[1]);
                break;
            case 3:
                out.coords[0] = e.e[0] - xdrop;
                out.coords[1] = e.l[1] - e.e[1];
                break;
            default:
                break;
        }
        return out;
    }
};

struct ReadOverlapExtraHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << e.dir << "\t" << e.sfx << "\t" << e.dirT << "\t" << e.sfxT << "\t" << static_cast<int>(e.transpose) << "\t" << static_cast<int>(e.rc) << "\t" << e.b[0] << "\t" << e.e[0] << "\t" << e.l[0] << "\t" << e.b[1] << "\t" << e.e[1] << "\t" << e.l[1];
    }
};

struct ReadOverlapCheckHandler
{
    template <typename c, typename t>
    void save(std::basic_ostream<c,t>& os, const ReadOverlap& e, int64_t row, int64_t col)
    {
        os << e.dir << "\t" << static_cast<int>(e.transpose) << "\t" << e.b[0] << "\t" << e.e[0] << "\t" << e.l[0] << "\t" << e.b[1] << "\t" << e.e[1] << "\t" << e.l[1];
    }

};

#endif
