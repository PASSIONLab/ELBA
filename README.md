# DISTAL: A Fast, Distributed Protein Sequence Aligner
=====

# Prerequisites
-----
1. Operating System.
  * DISTAL is tested and known to work on the following operating systems.
    *  SUSE Linux Enterprise Server 15.
    *  Ubuntu 14.10.
    *  MacOS.
    
2. GCC/G++ version 8.2.0 or above.

3. CMake 3.11 or above.


# Dependencies
-----
    
1. CombBLAS.
  * Download or clone CombBLAS from `https://bitbucket.org/berkeleylab/combinatorial-blas-2.0`
  * Export the path to this directory as an environment variable `COMBBLAS_HOME`
  * The following commands can be used to build and install CombBLAS
  ```
    cd $COMBBLAS_HOME/CombBLAS
    mkdir build
    mkdir install
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=../install ../
    make -j4
    make install         
  ```
3. SeqAn.
  * Download SeqAn `2.4.0` from `https://github.com/seqan/seqan/releases/tag/seqan-v2.4.0`
  * Extract this to a folder and create an environment variable, `SEQAN_HOME`, pointing to it. 
  * This is a header only library, so there's no need to build it.
  
  
# Build DISTAL
-----

To build DISTAL, you can clone or download the source from `