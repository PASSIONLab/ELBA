#ifndef SHARED_SEEDS_H_
#define SHARED_SEEDS_H_

#include "KmerOps.hpp"
#include <iostream>

using SeedPair = std::tuple<PosInRead, PosInRead>;

class Seed
{
public:
    Seed() : numshared(0) {}
    Seed(const Seed& rhs) : seeds(rhs.seeds), numshared(rhs.numshared) {}
    Seed(PosInRead begQ, PosInRead begT) : numshared(1)
    {
        std::get<0>(seeds[0]) = begQ;
        std::get<1>(seeds[0]) = begT;
    }

    int getnumstored() const { return std::min(MAX_SEEDS, numshared); }
    int getnumshared() const { return numshared; }

    void push(const SeedPair& seed)
    {
        if (numshared < MAX_SEEDS)
        {
            seeds[numshared] = seed;
        }

        numshared++;
    }

    const SeedPair* getseeds() const { return &seeds[0]; }

    Seed& operator=(const Seed& rhs)
    {
        Seed tmp(rhs);
        std::swap(seeds, tmp.seeds);
        numshared = tmp.numshared;
        return *this;
    }

    Seed& operator+=(const Seed& rhs)
    {
        push(rhs.seeds[0]);
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const Seed& o)
    {
        if (o.getnumstored() == 0)
        {
            os << "stored: {}";
        }
        else
        {
            os << "stored: {";
            int i;
            for (i = 0; i < o.getnumstored()-1; ++i)
                os << "(" << std::get<0>(o.seeds[i]) << ", " << std::get<1>(o.seeds[i]) << "), ";
            os << "(" << std::get<0>(o.seeds[i]) << ", " << std::get<1>(o.seeds[i]) << ")} :: numshared: " << o.getnumshared();
        }
        return os;
    }

    struct Semiring
    {
        static Seed id() { return Seed(); }
        static bool returnedSAID() { return false; }

        static Seed add(const Seed& lhs, const Seed& rhs)
        {
            Seed result(lhs);
            result += rhs;
            return result;
        }

        static Seed multiply(const PosInRead& lhs, const PosInRead& rhs)
        {
            return Seed(lhs, rhs);
        }

        static void axpy(PosInRead a, const PosInRead& x, Seed& y)
        {
            y = add(y, multiply(a, x));
        }
    };

    struct IOHandler
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t>& os, const Seed& o, uint64_t row, uint64_t col) { os << o; }
    };

private:
    SeedPair seeds[MAX_SEEDS];
    int numshared;
};

#endif
