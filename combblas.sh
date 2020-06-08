#!/bin/bash
git clone https://bitbucket.org/berkeleylab/combinatorial-blas-2.0
export COMBBLAS_HOME=combinatorial-blas-2.0/
cd $COMBBLAS_HOME/CombBLAS
mkdir build
mkdir install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
make -j4
make install 
