#ifndef SHARED_SEEDS_H_
#define SHARED_SEEDS_H_

#include "KmerOps.hpp"
#include <iostream>

using SeedPair = std::tuple<PosInRead, PosInRead>;

template <size_t MAXSEEDS>
class SharedSeeds
{
public:
    SharedSeeds() : seeds{}, numstored(0), numshared(0) {}
    SharedSeeds(const SharedSeeds& rhs) : seeds(rhs.seeds), numstored(rhs.numstored), numshared(rhs.numshared) {}
    SharedSeeds(PosInRead begQ, PosInRead begT) : seeds{}, numstored(1), numshared(1) { push(SeedPair(begQ, begT)); }
    /* SharedSeeds(PosInRead begQ, PosInRead begT) : seeds{}, numstored(1), numshared(1) { std::get<0>(seeds[0])} */

    int getnumstored() const { return numstored; }
    int getnumshared() const { return numshared; }

    void push(const SeedPair& seed)
    {
        numshared++;

        if (numstored < MAXSEEDS)
        {
            seeds[numstored++] = seed;
        }
    }

    const SeedPair* getseeds() const { return seeds.data(); }

    SharedSeeds& operator=(const SharedSeeds& rhs)
    {
        SharedSeeds tmp(rhs);
        std::swap(seeds, tmp.seeds);
        numshared = tmp.numshared;
        numstored = tmp.numstored;
        return *this;
    }

    SharedSeeds& operator+=(const SharedSeeds& rhs)
    {
        push(rhs.seeds.back());
        return *this;
    }

    friend std::ostream& operator<<(std::ostream& os, const SharedSeeds& o)
    {
        if (o.numstored == 0)
        {
            os << "stored: {}";
        }
        else
        {
            os << "stored: {";
            int i;
            for (i = 0; i < o.numstored-1; ++i)
                os << "(" << std::get<0>(o.seeds[i]) << ", " << std::get<1>(o.seeds[i]) << "), ";
            os << "(" << std::get<0>(o.seeds[i]) << ", " << std::get<1>(o.seeds[i]) << ")} :: numshared: " << o.numshared;
        }
        return os;
    }

    struct Semiring
    {
        static SharedSeeds id() { return SharedSeeds(); }
        static bool returnedSAID() { return false; }

        static SharedSeeds add(const SharedSeeds& lhs, const SharedSeeds& rhs)
        {
            SharedSeeds result(lhs);
            result += rhs;
            return result;
        }

        static SharedSeeds multiply(const PosInRead& lhs, const PosInRead& rhs)
        {
            return SharedSeeds(lhs, rhs);
        }

        static void axpy(PosInRead a, const PosInRead& x, SharedSeeds& y)
        {
            y = add(y, multiply(a, x));
        }
    };

    struct IOHandler
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t>& os, const SharedSeeds& o, uint64_t row, uint64_t col) { os << o; }
    };

private:
    std::array<SeedPair, MAXSEEDS> seeds;
    int numstored, numshared;
};

using Seed = SharedSeeds<2>;

#endif
