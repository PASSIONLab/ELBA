#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 1
#SBATCH -q debug
#SBATCH -J ecoli-gpuperl
#SBATCH --mail-user=gguidi@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 30:00
#SBATCH -A m2865_g

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export CURDIR=$PWD

#applications may performance better with --gpu-bind=none instead of --gpu-bind=single:1 
srun -n 64 -c 2 --cpu_bind=cores -G 1 --gpu-bind=none ${CURDIR}/build_release/./elba -i ${CURDIR}/ecoli_hifi_29x.fasta -k 31 --idxmap ecoli-idxmap -c 8605 --alph dna --af ecoli-cpu -s 1 -O 100000 --afreq 100000 --ca 15

export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export CURDIR=$PWD

srun -n 1 -c 128 --cpu_bind=cores -G 1 --gpu-bind=none ${CURDIR}/build_release/./elba -i ${CURDIR}/ecoli_hifi_29x.fasta -k 31 --idxmap ecoli-idxmap -c 8605 --alph dna --af ecoli-cpu -s 1 -O 100000 --afreq 100000 --ca 15
