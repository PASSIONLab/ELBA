# PASTIS-based diBELLA: Distributed Long Read to Long Read Alignment Using Sparse Matrices
=====

## Prerequisites

1. Operating System.
  * PASTIS is tested and known to work on the following operating systems.
    *  SUSE Linux Enterprise Server 15.
    *  Ubuntu 14.10.
    *  MacOS.
    
2. GCC/G++ version 8.2.0 or above.

3. CMake 3.11 or above.

4. Boost C++ Libraries.

## Dependencies
    
1. CombBLAS.
  * Download or clone CombBLAS from `https://bitbucket.org/berkeleylab/combinatorial-blas-2.0`.
  * Export the path to this directory as an environment variable `COMBBLAS_HOME`.
   ```
      export COMBBLAS_HOME=/path/to/combinatorial-blas-2.0
   ```
  * The following commands can be used to build and install CombBLAS:
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
  * Download SeqAn `2.4.0` from `https://github.com/seqan/seqan/releases/tag/seqan-v2.4.0` or clone the repository using the release tag as follows:
  ```
    git clone --branch seqan-v2.4.0 https://github.com/seqan/seqan.git
  ```
  * Create an environment variable, `SEQAN_HOME`, pointing to it:
  ```
    export SEQAN_HOME=/path/to/seqan
  ```
  * Checkout to `develop` branch of SeqAn `2.4.0`:
  ```
    cd $SEQAN_HOME
    git checkout develop
  ```
  * This is a header only library, so there's no need to build it.
  
## Use scripts to install dependencies
-----
  ```
    source modules.sh
    source combblas.sh
    cd ../../../
    source seqan.sh
    source utils.sh
    cd ../        
  ```
# Build diBELLA
-----

To build diBELLA, you can clone or download the source from here
  * mkdir build_release
  * cd build_release
  * cmake ..
  * cd ..
  * ./build.sh

# Run diBELLA
-----

You can run diBELLA in parallel by specifying the number of processes to the mpirun or mpiexec command. The number of processes must be perfect square value.

## Copy input data from m1982 (internal)
-----
  ```
    mkdir inputs
    lfs setstripe -c 72 -S 8M inputs/
    cd inputs/
    cp -r /project/projectdirs/m1982/bella-data/* .      
  ```
## From FASTQ to FASTA
-----
  ```
    cd ../seqtk
    ./seqtk seq -a <name>.fastq/fq > <name>.fa
  ```
## Ready to run
-----

The parameters and options of PASTIS are as follows:
- ```-i <string>```: Input FASTA file.
- ```-c <integer>```: Number of sequences in the FASTA file.
- ```--sc <integer>```: Seed count. ```[default: 2]```
- ```-k <integer>```: K-mer length.
- ```-s <integer>```: K-mers stride ```[default: 1]```
- ```--subs <integer>```: Number of substitute K-mers. 
- ```-g <integer>```: Gap open penalty (negative). ```[default: -11]```
- ```-e <integer>```: Gap extension penalty (negative). ```[default: -2]```
- ```-O <integer>```: Number of bytes to overlap when reading the input file in parallel. ```[default: 10000]```
- ```--na```: Do not perform alignment.
- ```--fa```: Smith-Waterman alignment.
- ```--xa <integer>```: X-drop alignment with the indicated drop value.
- ```--ba <integer>```: Banded alignment with the indicated band size.
- ```--af <string>```: Output file to write alignment information. 
- ```--idxmap <string>```: Output file for input sequences to ids used in PASTIS.
