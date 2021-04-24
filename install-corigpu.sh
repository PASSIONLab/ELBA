#!/bin/bash

source corigpu-env.sh
read lower upper delta

export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan
export LOGAN_HOME=$PWD/LoganGPU

echo ""
echo "Enviroment variables set and modules loaded."
echo ""

cd $LOGAN_HOME
mkdir build
cd build
cmake ..
make -j4

echo ""
echo "LOGAN installation completed."
echo ""

cd $COMBBLAS_HOME
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
cmake -DLOWER_KMER_FREQ=$lower -DUPPER_KMER_FREQ=$upper -DDELTACHERNOFF=$delta ..
make -j8

echo ""
echo "ELBA installation completed using upper [$upper], lower [$lower], and delta [$delta]."
echo ""
