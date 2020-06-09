# Internal Class Structure and Functionality
=====

# The Entry Point `main`
-----
The flow of diBELLA is as follows.
1. Initialize parallelism.
2. Parse arguments.
3. Create a unique job name (needed for logs)
4. Read local chunks from the Fasta file
5. Find what other processes needs to be contacted to send/recv sequences that may be necessary during the alignment phase.
6. Initiate async ISend/IRecv requests.
5. Discover k-mers in the local sequence set and create the `A` matrix.

# Fasta File Handling `ParallelFastaReader`
-----

