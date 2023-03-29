#!/bin/bash
#SBATCH -N 1
#SBATCH -C gpu
#SBATCH -G 4
#SBATCH -q regular
#SBATCH -J ecoli-cpu-gpu
#SBATCH --mail-user=gguidi@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH -t 45:00
#SBATCH -A m4341
#SBATCH --error=ecoli-elba-%j
#SBATCH --output=ecoli-elba-%j

#OpenMP settings:
export OMP_NUM_THREADS=1
export OMP_PLACES=threads
export OMP_PROC_BIND=true

export CURDIR=$PWD

#applications may performance better with --gpu-bind=none instead of --gpu-bind=single:1 
# CPU-based alignment
srun -n 1 -c 128 --cpu_bind=cores -G 1 --gpu-bind=none ${CURDIR}/build_release/./elba -i ${CURDIR}/ecoli_hifi_29x.fasta -k 31 --idxmap ecoli-idxmap -c 8605 --alph dna --af ecoli-cpu -s 1 -O 100000 --afreq 100000 --ca 15
# GPU-based alignment on 4 GPUs
srun -n 1 -c 128 --cpu_bind=cores -G 4 --gpu-bind=none ${CURDIR}/build_release/./elba -i ${CURDIR}/ecoli_hifi_29x.fasta -k 31 --idxmap ecoli-idxmap -c 8605 --alph dna --af ecoli-15-gpu-none -s 1 -O 100000 --afreq 100000 --ga 15

## this is 32 MPI processes per node, 1 thread per node, 4 GPUs

# #!/bin/bash
# #SBATCH -N 1
# #SBATCH -C gpu
# #SBATCH -G 4
# #SBATCH -q regular
# #SBATCH --mail-user=gguidi@cs.cornell.edu
# #SBATCH --mail-type=ALL
# #SBATCH -t 00:30:00
# #SBATCH -A m4341

# #OpenMP settings:
# export OMP_NUM_THREADS=1
# export OMP_PLACES=threads
# export OMP_PROC_BIND=spread

# #run the application:
# #applications may perform better with --gpu-bind=none instead of --gpu-bind=single:1 
# srun -n 32 -c 4 --cpu_bind=cores -G 4 --gpu-bind=single:1  myapp.x
