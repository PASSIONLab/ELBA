// Created by Saliya Ekanayake on 10/15/19.

#ifndef DISTAL_MATRIXENTRY_HPP
#define DISTAL_MATRIXENTRY_HPP

#include "../Types.hpp"
namespace distal{
  struct MatrixEntry{
    short cost;
    ushort offset;

    MatrixEntry():cost(0), offset(0){}
    MatrixEntry(short cost, ushort offset):cost(cost), offset(offset){}

    friend std::ostream& operator<<(std::ostream& os, const MatrixEntry& me){
      os << "(" << me.cost << ", " << me.offset << ")";
      return os;
    }

    MatrixEntry operator+(const MatrixEntry& me) const{
//      std::cout << "*** PLUS op was called";
      return me;
    }

    bool operator<(const MatrixEntry& me) const{
//      std::cout << "*** LESS THAN op was called";
      return cost < me.cost;
    }
  };
}
#endif //DISTAL_MATRIXENTRY_HPP
