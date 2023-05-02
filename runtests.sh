#!/bin/bash

diff -s test/ecoli.np1.txt <(mpirun -np 1 ./elba test/ecoli.fa)
diff -s gridrow1.txt <(seqkit seq -s test/ecoli.fa)
rm gridrow1.txt

diff -s test/ecoli.np4.txt <(mpirun -np 4 ./elba test/ecoli.fa)
diff -s gridrow1.txt <(seqkit seq -s test/ecoli.fa)
diff -s gridrow2.txt <(seqkit seq -s test/ecoli.fa)
rm gridrow{1..2}.txt

diff -s test/ecoli.np9.txt <(mpirun -np 9 ./elba test/ecoli.fa)
diff -s gridrow1.txt <(seqkit seq -s test/ecoli.fa)
diff -s gridrow2.txt <(seqkit seq -s test/ecoli.fa)
diff -s gridrow3.txt <(seqkit seq -s test/ecoli.fa)
rm gridrow{1..3}.txt

diff -s test/ecoli.np16.txt <(mpirun -np 16 ./elba test/ecoli.fa)
diff -s gridrow1.txt <(seqkit seq -s test/ecoli.fa)
diff -s gridrow2.txt <(seqkit seq -s test/ecoli.fa)
diff -s gridrow3.txt <(seqkit seq -s test/ecoli.fa)
diff -s gridrow4.txt <(seqkit seq -s test/ecoli.fa)
rm gridrow{1..4}.txt
