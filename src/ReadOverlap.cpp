#include "ReadOverlap.hpp"
#include <limits>
#include <assert.h>

static constexpr int MAX_INT = std::numeric_limits<int>::max();

int intplus(int a, int b)
{
    return (a == MAX_INT || b == MAX_INT)? MAX_INT : a + b;
}

void ReadOverlap::SetPathInf() { sfxpath[0] = sfxpath[1] = sfxpath[2] = sfxpath[3] = MAX_INT; }

ReadOverlap::ReadOverlap() : sfx(0), dir(-1), score(-1), count(1), transpose(false), passed(false) { SetPathInf(); }

ReadOverlap::ReadOverlap(int count) : ReadOverlap() { this->count = count; }

ReadOverlap::ReadOverlap(const ReadOverlap& rhs)
    : sfx(rhs.sfx), sfxT(rhs.sfxT), dir(rhs.dir), dirT(rhs.dirT), score(rhs.score), count(rhs.count),
      transpose(rhs.transpose), passed(rhs.passed), rc(rhs.rc)
{
    b[0] = rhs.b[0]; b[1] = rhs.b[1];
    e[0] = rhs.e[0]; e[1] = rhs.e[1];
    l[0] = rhs.l[0]; l[1] = rhs.l[1];

    for (int i = 0; i < 4; ++i)
        sfxpath[i] = rhs.sfxpath[i];

    begQs[0] = rhs.begQs[0];
    begQs[1] = rhs.begQs[1];

    begTs[0] = rhs.begTs[0];
    begTs[1] = rhs.begTs[1];
}

bool ReadOverlap::is_invalid() const { return (dir == -1); }

bool ReadOverlap::arrows(int& t, int& h) const
{
    if (is_invalid()) return false;

    t = (dir >> 1) & 1;
    h = dir & 1;

    return true;
}
