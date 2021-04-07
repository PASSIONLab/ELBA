// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "SR.hpp"
#include "CC.h"

/*! Namespace declarations */
using namespace combblas;

std::vector<std::string> 
CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu)
{
    // PSpMat<dibella::CommonKmers>::MPI_DCCols ST = S;
    // ST.Transpose();

    // if(!(ST == S)) tu.print_str("Error: matrix S isn't symmetric"); // The == overloading might make this wrong, fix it so it reflects reality

    uint nnz;
    PSpMat<dibella::CommonKmers>::MPI_DCCols S1 = S;

    /* Gonna iterate on the neighbor matrix is empty (no more neighbors to concatenate) */
    do
    {
        /* Find two-hops neighbors and concatenate sequences
         * N = S^2
         */

        // PSpMat<dibella::CommonKmers>::MPI_DCCols N = Mult_AnXBn_DoubleBuff<MinPlusSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(S, S1);
        PSpMat<dibella::CommonKmers>::MPI_DCCols N = Mult_AnXBn_DoubleBuff<MinPlusSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(S, S1);

        N.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);

        N.ParallelWriteMM("matrixN.mm", true, dibella::CkOutputMMHandler()); 
        tu.print_str("Matrix N = S*N: ");
        N.PrintInfo();
    #ifdef DIBELLA_DEBUG

    #endif

        // Do something to concantenate sequences based on S an N;

        // N is now the new matrix to be multiply S with and once is empty we can stop and print contigs?
        S1 = N; 

        nnz = N.getnnz(); 
            
    } while (nnz != 0);

    std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    S.ParallelWriteMM("matrixSc.mm", true, dibella::CkOutputMMHandler()); 
    tu.print_str("Matrix S: ");
    S.PrintInfo();

    return myContigSet;
}