#!/bin/bash

module unload PrgEnv-intel
module load PrgEnv-gnu

LOWER_KMER_FREQ=$1
UPPER_KMER_FREQ=$2
DELTACHERNOFF=$3
FUZZ=$4

module load cmake
module load boost

cd $HOME/software/mappers/diBELLA.2D

export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan
export MYPROJDIR=$PWD

rm -rf build_release
mkdir -p build_release
cd build_release

cmake .. -DLOWER_KMER_FREQ=$LOWER_KMER_FREQ -DUPPER_KMER_FREQ=$UPPER_KMER_FREQ -DDELTACHERNOFF=$DELTACHERNOFF -DFUZZ=$FUZZ
make -j8
cp dibella $BINARIES
cd $HOME
