#!/bin/bash

diff -s test/ecoli.np1.txt <(mpirun -np 1 ./elba test/ecoli.fa)
diff -s test/ecoli.np4.txt <(mpirun -np 4 ./elba test/ecoli.fa)
diff -s test/ecoli.np9.txt <(mpirun -np 9 ./elba test/ecoli.fa)
diff -s test/ecoli.np16.txt <(mpirun -np 16 ./elba test/ecoli.fa)
