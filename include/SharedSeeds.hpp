#ifndef SHARED_SEEDS_H_
#define SHARED_SEEDS_H_

#include "KmerOps.hpp"
#include <iostream>
#include <algorithm>

class SharedSeeds
{
public:
    SharedSeeds() : numshared(0) {}
    SharedSeeds(const SharedSeeds& rhs) : seeds(rhs.seeds), numshared(rhs.numshared) {}
    SharedSeeds(PosInRead begQ, PosInRead begT) : numshared(1)
    {
        std::get<0>(seeds[0]) = begQ;
        std::get<1>(seeds[0]) = begT;
    }

    int getnumstored() const { return std::min(2, numshared); }
    int getnumshared() const { return numshared; }

    const std::tuple<PosInRead, PosInRead>* getseeds() const { return &seeds[0]; }

    SharedSeeds& operator=(SharedSeeds rhs)
    {
        std::swap(seeds, rhs.seeds);
        numshared = rhs.numshared;
        return *this;
    }

    SharedSeeds& operator+=(const SharedSeeds& rhs)
    {
        if (numshared < 2)
        {
            seeds[numshared] = rhs.seeds[0];
        }

        numshared += rhs.numshared;

        return *this;
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
            SharedSeeds result(lhs, rhs);
            return result;
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

    struct IOHandlerBrief
    {
        template <typename c, typename t>
        void save(std::basic_ostream<c,t>& os, const SharedSeeds& o, uint64_t row, uint64_t col)
        {
            os << o.getnumstored() << "\t" << o.getnumshared();
        }
    };

    friend std::ostream& operator<<(std::ostream& os, const SharedSeeds& o)
    {
        int seedstoprint = std::min(static_cast<int>(o.getnumstored()), 2);

        os << "{";
        for (int i = 0; i < seedstoprint; ++i)
        {
            os << "(" << std::get<0>(o.seeds[i]) << "," << std::get<1>(o.seeds[i]) << "),";
        }

        os << o.getnumshared() << "}";

        return os;
    }

private:
    /*
     * @seeds is an array of size 2 whose elements (of type std::tuple<PosInRead, PosInRead>)
     * are ordered tuples of read positions.
     */
    std::tuple<PosInRead, PosInRead> seeds[2];
    int numshared;
};

std::unique_ptr<CT<SharedSeeds>::PSpParMat>
create_seed_matrix(CT<PosInRead>::PSpParMat& A, CT<PosInRead>::PSpParMat& AT);

#endif
