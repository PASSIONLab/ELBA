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
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt --idxmap idxmap.txt --subs
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt --idxmap idxmap.txt --subs
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap idxmap.txt --af align.txt --fa
#mpirun -np 1 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt


#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap idxmap.txt --sc 1 --na --of overlap.txt 2>&1 | tee na_out.txt
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap idxmap.txt --sc 1 --af xa_align.txt --xa 49 2>&1 | tee xa_out.txt
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap idxmap.txt --sc 1 --af ba_align.txt --ba 5 2>&1 | tee ba_out.txt
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap idxmap.txt --sc 1 --af fa_align.txt --fa --jp fa --lf 10 2>&1 | tee fa_out.txt
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/shuffled_100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap shuffled_idxmap.txt --sc 1 --af fa_shuffled_align.txt --fa --jp fa_shuffled --lf 10 2>&1 | tee fa_shuffled_out.txt

#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/1k/1000_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 1000 -O 1000 -k 6 -s 1 --idxmap idxmap.txt --sc 1 --af fa_align.txt --fa 2>&1 | tee fa_out.txt


# DEBUG runs
#mpirun -np 4 ./cmake-build-debug/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 6 -s 1 --of overlaps.txt --idxmap idxmap.txt
#fake run
#mpirun -np 1 ./cmake-build-debug/pisa -i /Users/esaliya/sali/data/scope/uniqs/10/fake_10_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 10 -O 1000 -k 1 -s 1 --of overlaps.txt --idxmap idxmap.txt

#mpirun -np 4 ./cmake-build-debug/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/shuffled_100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap shuffled_idxmap.txt --sc 1 --af fa_shuffled_align.txt --fa --jp fa_shuffled --lf 10 --subs 5 2>&1 | tee fa_shuffled_out.txt
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/shuffled_100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap shuffled_idxmap.txt --sc 1 --af fa_shuffled_align.txt --fa --jp fa_shuffled --lf 10 --subs 5 2>&1 | tee fa_shuffled_out.txt
#mpirun -np 4 ./cmake-build-release/pisa -i /Users/esaliya/sali/data/scope/uniqs/100/shuffled_100_of_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa -c 100 -O 1000 -k 6 -s 1 --idxmap shuffled_idxmap.txt --sc 1 --af fa_shuffled_align.txt --fa --jp fa_shuffled --lf 10 2>&1 | tee fa_shuffled_out.txt

#cog runs
#in_dir=/Users/esaliya/sali/data/cog/uniqs/shuffled
#in_file=74470_len_lte_100_in_shuffled_1769181_unique_of_1785722_prot2003-2014.fa
#mpirun -np 4 ./cmake-build-asan/pisa -i $in_dir/$in_file -c 74470 -O 10000 -k 6 -s 1 --idxmap cog_shuffled_idxmap.txt --sc 1  --subs 6 --af cog_fa_shuffled_align.txt --fa --jp cog_fa_shuffled --lf 10 2>&1 | tee cog_fa_shuffled_out.txt

#SEGFAULT when subk is increased to over 100 with SCOPe data
in_dir=/Users/esaliya/sali/data/scope/uniqs/all
in_file=10000_of_shuffled_77040_unique_of_243813_astral-scopedom-seqres-gd-all-2.07-stable.fa
mpirun -np 4 ./cmake-build-debug/pisa -i $in_dir/$in_file -c 10000 -O 10000 -k 6 -s 1 --idxmap shuffled_idxmap.txt --sc 1 --na --of na_overlap.txt --jp na_shuffled --subs 500 --lf 10000 2>&1 | tee na_shuffled_out.txt