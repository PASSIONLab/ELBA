# diBELLA: Distributed Long Read to Long Read Alignment Using Sparse Matrices
=====

## Prerequisites

1. Operating System.
  * diBELLA is tested and known to work on the following operating systems.
    *  SUSE Linux Enterprise Server 15.
    *  Ubuntu 14.10.
    *  MacOS.
    
2. GCC/G++ version 8.2.0 or above.

3. CMake 3.11 or above.

## Dependencies
    
1. CombBLAS.
  * Download or clone CombBLAS from `https://github.com/PASSIONLab/CombBLAS.git`.
  * Export the path to this directory as an environment variable `COMBBLAS_HOME`.
   ```
      git clone https://github.com/PASSIONLab/CombBLAS.git
      export COMBBLAS_HOME=$PWD
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
3. SeqAn (included in this repository).
  * Create an environment variable, `SEQAN_HOME`, pointing to it:
  ```
    export SEQAN_HOME=/path/to/seqan
  ```
  * This is a header only library, so there's no need to build it.

# Build diBELLA
-----

To build diBELLA, you can use the following commands:
  ```
    mkdir build_release
    cd build_release
    cmake -DCMAKE_C_COMPILER=/usr/local/bin/gcc-10 -DCMAKE_CXX_COMPILER=/usr/local/bin/g++-10 ..
    make -j4  
  ```
Default macro definition in CMakeFiles.txt:
  ```
    #define MAX_KMER_SIZE 32
    #define ERR_THRESHOLD 2
    #define MAX_NUM_READS 8
  ```
Based on the dataset, one might want to change the above definitions. *MAX_NUM_READS*: reliable k-mer upper bound (8 works for E. coli (Sample) 30X and 4 for Human CCS), *ERR_THRESHOLD*: reliable k-mer lower bound.

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

The parameters and options of diBELLA are as follows:
- ```-i <string>```: Input FASTA file.
- ```-c <integer>```: Number of sequences in the FASTA file.
- ```--sc <integer>```: Seed count. ```[default: 2]```
- ```-k <integer>```: K-mer length.
- ```-s <integer>```: K-mers stride. ```[default: 1]```
- ```--ma <integer>```: Base match score (positive). ```[default: 1]```
- ```--mi <integer>```: Base mismatch score (negative). ```[default: -1]```
- ```-g <integer>```: Gap open penalty (negative). ```[default: 0]```
- ```-e <integer>```: Gap extension penalty (negative). ```[default: -1]```
- ```-O <integer>```: Number of bytes to overlap when reading the input file in parallel. ```[default: 10000]```
- ```--afreq <integer>```: Alignment write frequency.
- ```--na```: Do not perform alignment.
- ```--fa```: Full Smith-Waterman alignment.
- ```--xa <integer>```: X-drop alignment with the indicated drop value.
- ```--ba <integer>```: Banded alignment with the indicated band size.
- ```--of <string>```: Overlap file.
- ```--af <string>```: Output file to write alignment information. 
- ```--idxmap <string>```: Output file for input sequences to ids used in diBELLA.
- ```--alph <dna|protein>```: Alphabet.
