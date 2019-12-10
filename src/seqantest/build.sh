#!/bin/bash
rm -rf a.mavx2.out
g++-9 -mavx2 -I/Users/esaliya/sali/software/seqan-v2.4.0/include -std=c++14 seqanvec.cpp -o a.mavx2.out
