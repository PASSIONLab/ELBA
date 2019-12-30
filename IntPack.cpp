// Created by Saliya Ekanayake on 10/17/19.
#include <iostream>

int64_t Pack(short a, unsigned short b)
{
  return (int64_t)((((uint64_t)a)<<16)+(uint32_t)b);
}

short UnpackA(int x)
{
  return (short)(((uint64_t)x)>>16);
}

unsigned short UnpackB(int x)
{
  return (unsigned short)(((uint64_t)x)&0xffff);
}
int main(int argc, char** argv){
  short x = -210;
  unsigned short y = 2450;
  int z = Pack(x, y);
  short xp = UnpackA(z);
  unsigned short yp = UnpackB(z);
  std::cout << xp << " " << yp << std::endl;
}