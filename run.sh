#!/bin/bash
#mpirun -np 4 ./cmake-build-release/lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5_test.fasta -c 5 -O 1000 -k 6 -s 1
#mpirun -np 4 ./cmake-build-debug/lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5_test.fasta -c 5 -O 1000 -k 6 -s 1
#mpirun -np 9 ./cmake-build-debug/lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 6 -s 1
#mpirun -np 9 ./cmake-build-release/lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 6 -s 1

mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 1 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5.fasta -c 5 -O 1000 -k 2 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_6.fasta -c 6 -O 1000 -k 2 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_7.fasta -c 7 -O 1000 -k 2 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_tmp.txt -c 5 -O 1000 -k 2 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_8.fasta -c 8 -O 1000 -k 2 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_10.fasta -c 10 -O 1000 -k 2 -s 1
#mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_50.fasta -c 50 -O 1000 -k 2 -s 1
#mpirun -np 16 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 2 -s 1

#mpirun -np 16 ./cmake-build-release/lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 6 -s 1


# Astral (SCOPe) sequences
#mpirun -np 4 ./cmake-build-release/lbl_dal -i /Users/esaliya/sali/data/scope/5000_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 5000 -O 1000 -k 6 -s 1
