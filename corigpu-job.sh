#!/bin/bash
#SBATCH -A m2865
#SBATCH -C gpu
#SBATCH -N 1
#SBATCH -q regular
#SBATCH -t 30:00
#SBATCH --ntasks-per-node=1
#SBATCH -c 40
#SBATCH --gpus-per-node=8

# Despite what its name suggests, --gpus-per-task in the examples below only counts the number of GPUs to allocate to the job; it does not enforce any binding or affinity of GPUs to CPUs or tasks.

export OMP_NUM_THREAD=40
export SLURM_CPU_BIND="cores"

srun -n 1 ~/ELBA/build_release/./elba -i ~/ELBA/sub_sample_ecoli_hifi.fasta -k 31 --idxmap ecoli.hifi.idxmap -c 1000 --alph dna --af ecoli.hifi.result -s 1 -O 100000 --afreq 100000 --ga 15

rm readNameMap_* elba_rank_* contiglog_*

srun -n 1 ~/ELBA/build_release/./elba -i ~/ELBA/sub_sample_ecoli_hifi.fasta -k 31 --idxmap ecoli.hifi.idxmap -c 1000 --alph dna --af ecoli.hifi.result -s 1 -O 100000 --afreq 100000 --ca 15

rm readNameMap_* elba_rank_* contiglog_*

