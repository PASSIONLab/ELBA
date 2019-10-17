// Created by Saliya Ekanayake on 10/16/19.

struct ME{
  short cost;
  unsigned short offset;

  MatrixEntry():cost(0), offset(0){}
  MatrixEntry(short cost, ushort offset):cost(cost), offset(offset){}

  friend std::ostream& operator<<(std::ostream& os, const MatrixEntry& me){
    os << "(" << me.cost << ", " << me.offset << ")";
    return os;
  }

  MatrixEntry operator+(const MatrixEntry& me) const{
    MatrixEntry me2(cost,offset);
    return me2;
  }

  bool operator<(const MatrixEntry& me) const{
    return true;
  }
};
int main(int argc, char** argv){

}