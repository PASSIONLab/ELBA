#!/bin/bash
module swap PrgEnv-intel PrgEnv-gnu && module load cmake && module load boost && module load python
export COMBBLAS_HOME=/global/cscratch1/sd/gguidi/IPDPS2021/diBELLA/CombBLAS
export BLOOM_HOME=/global/cscratch1/sd/gguidi/IPDPS2021/diBELLA/src/libbloom/
export SEQAN_HOME=/global/cscratch1/sd/gguidi/IPDPS2021/diBELLA/seqan
