// Created by Saliya Ekanayake on 10/15/19.

#ifndef LBL_PISA_MATRIXENTRY_HPP
#define LBL_PISA_MATRIXENTRY_HPP

#include "../Types.hpp"
namespace pisa{
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
      std::cout << "*** PLUS op was called";
      MatrixEntry me2(cost,offset);
      return me2;
    }

    bool operator<(const MatrixEntry& me) const{
      std::cout << "*** LESS THAN op was called";
      return true;
    }
  };
}
#endif //LBL_PISA_MATRIXENTRY_HPP
