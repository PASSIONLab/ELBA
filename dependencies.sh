#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu
module load openmpi
module load cmake/3.14.4
git clone https://bitbucket.org/berkeleylab/combinatorial-blas-2.0
export COMBBLAS_HOME=combinatorial-blas-2.0/
cd $COMBBLAS_HOME/CombBLAS
mkdir build
mkdir install
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
make -j4
make install
cd ../../../
wget https://github.com/seqan/seqan/releases/download/seqan-v2.4.0/seqan-library-2.4.0.zip
unzip seqan-library-2.4.0
rm seqan-library-2.4.0.zip
export SEQAN_HOME=seqan-library-2.4.0/
