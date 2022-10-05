#!/bin/bash
#SBATCH -A m2865
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -t 5:00
#SBATCH --ntasks-per-node=16
#SBATCH -c 1
#SBATCH --gpus-per-node=8

# Despite what its name suggests, --gpus-per-task in the examples below only counts the number of GPUs to allocate to the job; it does not enforce any binding or affinity of GPUs to CPUs or tasks.

export OMP_NUM_THREAD=40
export SLURM_CPU_BIND="cores"
srun -n 16 /global/cscratch1/sd/gguidi/ELBA-GPU/ELBA/build_release/./elba -i /global/cscratch1/sd/gguidi/ELBA-GPU/ELBA/ecoli_hifi_29x.fasta -k 31 --idxmap ecoli.hifi.idxmap -c 8605 --alph dna --af ecoli.hifi.result -s 1 -O 100000 --afreq 100000 --xa 15
