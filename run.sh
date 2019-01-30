#!/bin/bash
mpirun -np 4 ./lbl_dal -i /Users/esaliya/sali/data/metaclust_50/metaclust_50_head_5.fasta -c 10 -O 1000 -k 6 -s 1
