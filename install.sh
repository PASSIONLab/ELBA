#!/bin/bash
export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan

git clone https://github.com/PASSIONLab/CombBLAS.git
cd $COMBBLAS_HOME/CombBLAS
mkdir build
mkdir install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
make -j8
make install

echo ""
echo "CombBLAS installation completed."
echo ""

cd $COMBBLAS_HOME
mkdir build_release
cd build_release
cmake ..
make -j8

echo ""
echo "ELBA installation completed."
echo ""
