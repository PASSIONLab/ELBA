#!/bin/bash

module swap PrgEnv-intel PrgEnv-gnu && module load cmake && module load boost && module load python
export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan

git clone https://github.com/gabe-raulet/CombBLAS
cd $COMBBLAS_HOME/CombBLAS
git checkout from-scratch
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

echo "diBELLA 2D installation completed."
echo ""
