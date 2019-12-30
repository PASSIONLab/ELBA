#!/bin/bash
#srun -n 4 valgrind ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_5_test.fasta -c 5 -O 1000 -k 6 -s 1 2>&1 | tee out.txt
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_5_test.fasta -c 5 -O 1000 -k 6 -s 1


#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_5.fasta -c 5 -O 1000 -k 2 -s 1
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_6.fasta -c 6 -O 1000 -k 2 -s 1
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_7.fasta -c 7 -O 1000 -k 2 -s 1
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_tmp.txt -c 5 -O 1000 -k 2 -s 1
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_8.fasta -c 8 -O 1000 -k 2 -s 1
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_10.fasta -c 10 -O 1000 -k 2 -s 1
#srun -n 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_50.fasta -c 50 -O 1000 -k 2 -s 1
#srun -n 16 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 2 -s 1


#srun -n 9 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_100.fasta -c 100 -O 1000 -k 1 -s 1 2>&1 | tee out.txt

#srun -n 25 valgrind ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 6 -s 1 2>&1 | tee out.txt
#srun -n 16 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_4000.fasta -c 4000 -O 1000 -k 6 -s 1



#srun -n 16 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_25000.fasta -c 25000 -O 1000 -k 4 -s 1

#srun -n 25 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_999991.fasta -c 999991 -O 1000 -k 6 -s 1

#PERF RUNS
#srun -n 16 -c 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_25000.fasta -c 25000 -O 1000 -k 6 -s 1 2>&1 | tee m50.25k.k6.s1.1x16x1.c4.txt
#srun -n 16 -c 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_100k.fasta -c 100000 -O 1000 -k 6 -s 1 2>&1 | tee m50.100k.k6.s1.1x16x1.c4.txt
#srun -n 16 -c 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_250k.fasta -c 100000 -O 1000 -k 6 -s 1 2>&1 | tee m50.250k.k6.s1.1x16x1.c4.txt
#srun -n 16 -c 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_500k.fasta -c 100000 -O 1000 -k 6 -s 1 2>&1 | tee m50.500k.k6.s1.1x16x1.c4.txt
srun -n 16 -c 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_1M.fasta -c 1000000 -O 1000 -k 6 -s 1 2>&1 | tee m50.1M.k6.s1.1x16x1.c4.txt
srun -n 16 -c 4 ./lbl_dal -i /global/homes/e/esaliya/sali/data/metaclust_50/metaclust_50_head_5M.fasta -c 5000000 -O 1000 -k 6 -s 1 2>&1 | tee m50.5M.k6.s1.1x16x1.c4.txt
