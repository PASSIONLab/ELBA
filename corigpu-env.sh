module purge
module load cgpu
module load pgi
module load cmake/3.14.4
module load cuda 
module load openmpi
module load boost
module load python

export COMBBLAS_HOME=$PWD
export BLOOM_HOME=$PWD/src/libbloom/
export SEQAN_HOME=$PWD/seqan
export MYPROJDIR=$PWD
