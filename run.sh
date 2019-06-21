#!/bin/bash
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5_test.fasta -c 5 -O 1000 -k 6 -s 1
#mpirun -np 4 ./cmake-build-debug/pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5_test.fasta -c 5 -O 1000 -k 6 -s 1
#mpirun -np 9 ./cmake-build-debug/pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 6 -s 1
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 6 -s 1

#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 1 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5.fasta -c 5 -O 1000 -k 2 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_6.fasta -c 6 -O 1000 -k 2 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_7.fasta -c 7 -O 1000 -k 2 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_tmp.txt -c 5 -O 1000 -k 2 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_8.fasta -c 8 -O 1000 -k 2 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_10.fasta -c 10 -O 1000 -k 2 -s 1
#mpirun -np 4 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_50.fasta -c 50 -O 1000 -k 2 -s 1
#mpirun -np 16 ./pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 2 -s 1

#mpirun -np 16 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 6 -s 1


# Astral (SCOPe) sequences
mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt 
#mpirun -np 1 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt


# DEBUG runs
#mpirun -np 4 ./cmake-build-debug/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt
