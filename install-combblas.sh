#!/bin/bash

git clone https://github.com/PASSIONLab/CombBLAS
cd $COMBBLAS_HOME/CombBLAS
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
make -j8
make install
