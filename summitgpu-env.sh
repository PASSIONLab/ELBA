module load gcc/9.1.0
# It should be cmake/3.18
module load cmake
module load cuda
module load boost
module load python

# Otherwise it uses gcc-4.8
export CXX=g++
export CCC=gcc
