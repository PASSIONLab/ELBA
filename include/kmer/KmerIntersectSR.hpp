// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef ELBA_KMERINTERSECTSR_HPP
#define ELBA_KMERINTERSECTSR_HPP

#include "../ParallelOps.hpp"
#include "../Defines.hpp"

#include <cstdlib>
#include <cmath>
#include <functional>

std::pair<int, int> distance(const std::pair<PosInRead, PosInRead>& l, const std::pair<PosInRead, PosInRead>& r) {
	return {std::abs(static_cast<int>(l.first - r.first)), std::abs(static_cast<int>(l.second - r.second))};
}

bool operator>(const std::pair<int, int>& l, const int& c) {
	if(l.first > c && l.second > c) return true;
	else return false;
}

namespace elba {
  template<typename IN, typename OUT>
  struct KmerIntersect {
    static OUT id() {
      OUT a;
      return a;
    }

    static bool returnedSAID() { return false; }

    static OUT add(const OUT &arg1, const OUT &arg2)
    {
  #ifdef TWOSEED
      OUT res(arg1.count + arg2.count);

      res.first.first   = arg1.first.first; // pos of k1 on read1
      res.second.first  = arg2.first.first; // pos of k2 on read1

      res.first.second  = arg1.first.second; // pos of k1 on read2
      res.second.second = arg2.first.second; // pos of k2 on read2

      return res;
  #else
      OUT res(arg1.count);
      res.pos = arg1.pos;

      std::vector<std::pair<PosInRead, PosInRead>> kmertobeinserted; //(arg1->pos.size()); GGGG: techinically i don't need size

      for(int i = 0; i < arg2.pos.size(); ++i)	
        for(int j = 0; j < arg1.pos.size(); ++j)
        {
          auto kmer1 = arg1.pos[j];
          auto kmer2 = arg2.pos[i];

          if(distance(kmer1, kmer2) > KLEN)
            kmertobeinserted.push_back(kmer2);
        }

      for (int i = 0; i < kmertobeinserted.size(); i++)
      {
        res.count	+= kmertobeinserted.size();
        res.pos.insert(res.pos.end(), kmertobeinserted.begin(), kmertobeinserted.end());
      } 

      return res;
  #endif
    }

    static OUT multiply(const IN &arg1, const IN &arg2)
    {
      OUT a;

  #ifdef TWOSEED
      a.first.first  = arg1; // pos of k1 on read1
      a.first.second = arg2; // pos of k1 on read2
  #else
      std::pair<IN, IN> mypair{std::make_pair(arg1, arg2)};
      a.pos.push_back(mypair);
  #endif
  
      return a;
    }

    static void axpy(IN a, const IN &x, OUT &y) {
      y = add(y, multiply(a, x));
    }

    static MPI_Op mpi_op() {
      static MPI_Op mpiop;
      static bool exists = false;
      if (exists)
        return mpiop;
      else {
        MPI_Op_create(MPI_func, true, &mpiop);
        exists = true;
        return mpiop;
      }
    }

    static void
    MPI_func(void *invec, void *inoutvec, int *len, MPI_Datatype *datatype) {
      for (int i = 0; i < *len; ++i) {
        *((OUT) inoutvec + i) = add(*((OUT) invec + i), *((OUT) inoutvec + 1));
      }

    }
  };
}
#endif //ELBA_KMERINTERSECTSR_HPP
