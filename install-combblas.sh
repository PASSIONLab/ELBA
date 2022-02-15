#!/bin/bash

git clone https://github.com/gabe-raulet/CombBLAS
cd $COMBBLAS_HOME/CombBLAS
git checkout dibella-branch
mkdir -p build
mkdir -p install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
make -j8
make install
