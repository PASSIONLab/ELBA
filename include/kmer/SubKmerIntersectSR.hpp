// Created by Saliya Ekanayake on 10/15/19.

#ifndef DISTAL_SUBKMERINTERSECTSR_HPP
#define DISTAL_SUBKMERINTERSECTSR_HPP

#include "../ParallelOps.hpp"
#include "MatrixEntry.hpp"

#include "../ParallelOps.hpp"
namespace distal {
  template<typename IN, typename OUT>
  struct SubKmerIntersect {
    static OUT id() {
      OUT a;
      return a;
    }

    static bool returnedSAID() { return false; }

    static OUT add(const OUT &arg1, const OUT &arg2) {
      if (arg1.cost < arg2.cost) {
        OUT ret(arg1.cost, arg1.offset);
        return ret;
      } else {
        OUT ret(arg2.cost, arg2.offset);
        return ret;
      }
    }

    static OUT multiply(const IN &arg1, const IN &arg2) {
      OUT a(arg2.cost, arg1.offset);
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
#endif //DISTAL_SUBKMERINTERSECTSR_HPP
