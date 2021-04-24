#!/bin/bash
#SBATCH -A m2865
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 1:00:00
#SBATCH --ntasks-per-node=32
#SBATCH -c 2
#SBATCH --gpus-per-node=8

# Despite what its name suggests, --gpus-per-task in the examples below only counts the number of GPUs to allocate to the job; it does not enforce any binding or affinity of GPUs to CPUs or tasks.

export SLURM_CPU_BIND="cores"
srun -n 32 /global/cscratch1/sd/gguidi/diBELLA.2D/./dibella -i /global/cscratch1/sd/gguidi/diBELLA.2D/ecoli.hifi.29x.fasta -k 31 --idxmap ecoli.hifi.29x.idxmap -c 95514 --alph dna --af ecoli.hifi.29x.result -s 1 -O 100000 --afreq 100000 --xa 15
