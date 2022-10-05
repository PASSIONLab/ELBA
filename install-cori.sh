#!/bin/bash

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
mkdir build_release
cd build_release
cmake -DLOWER_KMER_FREQ=20 -DUPPER_KMER_FREQ=30 -DDELTACHERNOFF=0.7 ..
make -j8

echo ""
echo "ELBA installation completed using upper [30], lower [20], and delta [0.7]."
echo ""

