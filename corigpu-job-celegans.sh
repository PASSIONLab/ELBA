#!/bin/bash
#SBATCH -A m2865
#SBATCH -N 16
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=16
#SBATCH -c 2
#SBATCH --gpus-per-node=8

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export SLURM_CPU_BIND="cores"
srun -n 16 /global/cscratch1/sd/gguidi/sc-gpu/diBELLA.2D/build_release/./dibella -i /global/cscratch1/sd/gguidi/sc-gpu/diBELLA.2D/celegans40x.fa -k 17 --idxmap elba.celegans40x.skylake.nvidia.16.idxmap -c 4421593 --alph dna --af elba.celegans40x.skylake.nvidia.16.result -s 1 -O 100000 --afreq 100000 --xa 15
