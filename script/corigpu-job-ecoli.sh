#!/bin/bash
#SBATCH -A m2865
#SBATCH -C gpu
#SBATCH -q regular
#SBATCH -t 10:00
#SBATCH -N 1
#SBATCH --job-name=elba.logan.ecoli.ccs.29x.16
#SBATCH --ntasks-per-node=16
#SBATCH -c 2
#SBATCH --mail-user=gguidi@berkeley.edu
#SBATCH --mail-type=ALL
#SBATCH --gpus-per-node=8

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

# Despite what its name suggests, --gpus-per-task in the examples below only counts the number of GPUs to allocate to the job; it does not enforce any binding or affinity of GPUs to CPUs or tasks.

export SLURM_CPU_BIND="cores"
srun -n 16 /global/cscratch1/sd/gguidi/diBELLA.2D/build_release/./dibella -i /global/cscratch1/sd/gguidi/diBELLA.2D/ecoli_hifi_29x.fasta -k 31 --idxmap ecoli.hifi.29x.idxmap -c 17210 --alph dna --af ecoli.hifi.29x.result -s 1 -O 100000 --afreq 100000 --xa 15
