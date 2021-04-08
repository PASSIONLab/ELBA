// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "ContigEntry.hpp"
#include "Utils.hpp"
#include "SR.hpp"
#include "CC.h"

/*! Namespace declarations */
using namespace combblas;
typedef ContigSRing <dibella::ContigEntry, dibella::ContigEntry, dibella::ContigEntry> ContigSRing_t;

std::vector<std::string> 
CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu)
{    
    dibella::ContigEntry BNullVal;
    dibella::CommonKmers ANullVal;
    PSpMat<dibella::ContigEntry>::MPI_DCCols T;

    // T has dir, len of suffix and lenh, lenv of sequences -- this information is enough to chop the suffix sequence
    T = EWiseApply<dibella::ContigEntry, PSpMat<dibella::ContigEntry>::DCCols>(S, T, ContigEntrySR<dibella::CommonKmers, dibella::ContigEntry>(), 
                            ContigEntrySRP<dibella::CommonKmers, dibella::ContigEntry>(), false, true, ANullVal, BNullVal, false);

    // TODO: include sequences.

    PSpMat<dibella::ContigEntry>::MPI_DCCols nT = T; // neighbor T matrix

    /* Gonna iterate on the neighbor matrix is empty (no more neighbors to concatenate) */
    uint nnz;
    do
    {
        nT = Mult_AnXBn_DoubleBuff<ContigSRing_t, dibella::ContigEntry, PSpMat<dibella::ContigEntry>::DCCols>(T, nT);

        // nT.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);

        nT.ParallelWriteMM("matrixnT.mm", true, dibella::CeOutputMMHandler()); 
        tu.print_str("nT matrix = T*nT: ");
        nT.PrintInfo();
    #ifdef DIBELLA_DEBUG

    #endif

        nnz = nT.getnnz(); 
            
    } while (nnz != 0);

    std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    return myContigSet;
}