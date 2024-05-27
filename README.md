# ELBA
## Parallel String Graph Construction, Transitive Reduction, and Contig Generation for De Novo Genome Assembly

## Prerequisites

1. Operating System.
  * ELBA is tested and known to work on the following operating systems.
    *  SUSE Linux Enterprise Server 15.
    *  Ubuntu 14.10.
    *  MacOS.
    
2. GCC/G++ version 8.2.0 or above.

3. CMake 3.11 or above.

## Dependencies
    
0. STILL NEED TO EDIT THIS, WILL BE INCONSISTENT WITH ACTUAL REPO FOR A BIT.

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
    export BLOOM_HOME=src/libbloom/
  ```
  * This is a header only library, so there's no need to build it.

# Build ELBA
To build ELBA, you can use the following commands:
  ```
    mkdir build_release
    cd build_release
    cmake ..
    make -j4  
  ```
Default macro definition in CMakeFiles.txt:
  ```
    #define MAX_KMER_SIZE  32
    #define LOWER_KMER_FREQ 2
    #define UPPER_KMER_FREQ 8
  ```
Based on the dataset, one might want to change the above definitions. **UPPER_KMER_FREQ**: reliable k-mer upper bound (8 works for E. coli (Sample) 30X and 4 for Human 10X and C. elegans 40X that you can find [here](https://portal.nersc.gov/project/m1982/dibella.2d/inputs/)), **LOWER_KMER_FREQ**: reliable k-mer lower bound.

You can change the defaul setting at compile time when building using the following command instead of ```cmake ..```:
```
cmake -DLOWER_KMER_FREQ=<new-lower-bound> -DUPPER_KMER_FREQ=<new-upper-bound> .. 
```

# Run ELBA

You can run ELBA in parallel by specifying the number of processes to the mpirun or mpiexec command. The number of processes must be perfect square value.

## Input data samples
A few input data sets can be downloaded [here](https://portal.nersc.gov/project/m1982/dibella.2d/inputs/). If you have your own FASTQs, you can convert them into FASTAs using [seqtk](https://github.com/lh3/seqtk):

  ```
    cd ../seqtk
    ./seqtk seq -a <name>.fastq/fq > <name>.fa
  ```
A tiny example `ecsample-sub1.fa` can be found in this repository.

## Ready to run
The parameters and options of ELBA are as follows:
- ```-i <string>```: Input FASTA file.
- ```-x <integer>```: x-drop alignment threshold.
- ``` -A <integer>```: matching score.
- ```-B <integer>```: mismatch penalty.
- ```-G <integer>```:  gap penalty.
- ```-c <float> ```: bad read alignment cutoff.
- ```-o string```L output file name prefix "elba."
- ```-h```: help message.

## Run test program
You can run the test dataset ```ecsample-sub1.fa``` as follows on one node (it's too small to run on multiple nodes), this command runs ELBA using x-drop alignment and ```x = 5```:
```
export OMP_NUM_THREADS=1
mpirun -np 1 ./elba -i /path/to/ecsample-sub1.fa -k 17 --idxmap elba-test -c 135 --alph dna --of overlap-test --af alignment-test -s 1 -O 100000 --afreq 100000 --xa 5
```
To run on multiple nodes, for example on 4 nodes using 4 MPI rank/node, please download ```ecsample30x.fa``` from [here](https://portal.nersc.gov/project/m1982/dibella.2d/inputs/) and run as follows:
```
export OMP_NUM_THREADS=1
mpirun -np 16 ./elba -i /path/to/ecsample30x.fa -k 17 --idxmap elba-ecsample -c 16890 --alph dna --of overlap-ecsample --af alignment-ecsample -s 1 -O 100000 --afreq 100000 --xa 5
```
You need to use a perfect square number of processes to match our 2D decomposition. Recall ```-c``` should match the number of sequences in the input FASTA.

# Citation
To cite our work or to know more about our methods, please refer to:

> Giulia Guidi, Oguz Selvitopi, Marquita Ellis, Leonid Oliker, Katherine Yelick, Aydın Buluç. [Parallel String Graph Construction and Transitive Reduction for De Novo Genome Assembly](https://arxiv.org/pdf/2010.10055.pdf). Proceedings of the IPDPS, 2021.

> Giulia Guidi, Gabriel Raulet, Daniel Rokhsar, Leonid Oliker, Katherine Yelick, Aydın Buluç. [Distributed-Memory Parallel Contig Generation for De Novo
Long-Read Genome Assembly](https://arxiv.org/pdf/2207.04350.pdf). Proceedings of the ICPP, 2022.

Further design choices and results in terms of accuracy can be found here:

> Giulia Guidi, Marquita Ellis, Daniel Rokhsar, Katherine Yelick, Aydın Buluç. [BELLA: Berkeley Efficient Long-Read to Long-Read Aligner and Overlapper](https://drive.google.com/file/d/132i0RAKyIIWk_BEl1jpf9R_V5eVkKkxT/view). bioRxiv 464420; doi: https://doi.org/10.1101/464420. Proceedings of the SIAM ACDA, 2021.

# Copyright

diBELLA 2D: Parallel String Graph Construction and Transitive Reduction for De Novo Assembly (diBELLA 2D) Copyright (c) 2021, The Regents of the University of California, through Lawrence Berkeley National Laboratory (subject to receipt of any required approvals from the U.S. Dept. of Energy).  All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual Property Office at IPO@lbl.gov.

NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare derivative works, and perform publicly and display publicly, and to permit others to do so.
