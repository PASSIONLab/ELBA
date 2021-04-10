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
typedef ContigSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ContigSRing_t;

std::vector<std::string> 
CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu)
{    

    float balance = S.LoadImbalance();
    int64_t nnz   = S.getnnz();

    std::ostringstream outs;
    outs.str("");
    outs.clear();
    outs << "CreateContig::LoadBalance: " << balance << endl;
    outs << "CreateContig::nonzeros: "    << nnz     << endl;
    SpParHelper::Print(outs.str());

    int64_t nCC = 0;
    FullyDistVec<int64_t, int64_t> myLabelCC = CC(S, nCC);

	std::string myccoutput = myoutput + ".cc";
	myLabelCC.ParallelWrite(myccoutput, 1);

    uint64_t nContig = myLabelCC.Reduce(maximum<int64_t>(), (int64_t)0);
    nContig++; // because of zero based indexing for cluster

    std::stringstream ncc;
    ncc << "nContig: " << nContig << endl;
    SpParHelper::Print(ncc.str());

    std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    return myContigSet;
}

// PSpMat<dibella::CommonKmers>::MPI_DCCols T    = S; // (T) traversal matrix
// PSpMat<dibella::CommonKmers>::MPI_DCCols APSP = S; // (APSP) all-pairs shortest paths matrix
// dibella::CommonKmers defaultBVal; 
// do
// {   
//     PSpMat<dibella::CommonKmers>::MPI_DCCols P = // (P) power matrix
//         Mult_AnXBn_DoubleBuff<ContigSRing_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(S, T);

//     P.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);

//     APSP += P;

//     // update matrix
//     T = APSP;
//     nnz = P.getnnz();
            
// } while (nnz != 0);

// fname = "matrixAPSP-t-" + std::to_string(APSP.getnnz()) + ".mm";
// APSP.ParallelWriteMM(fname, true, dibella::CkOutputMMHandler());

// // I find the previous path
// PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(APSP, T, 
//     EqualBinaryOp<dibella::CommonKmers, dibella::CommonKmers>(), false, defaultBVal);
// I.Prune(ZeroUnaryOp<bool>(), true);
// // I remove the previous path so that APSP only has the max-hop one (contig)
// APSP = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(APSP, I, 
//     EWiseMulOp<dibella::CommonKmers, bool>(), true, true);
// // Prune zero-valued overhang
// APSP.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);
// fname = "matrixAPSP" + std::to_string(APSP.getnnz()) + ".mm";
// APSP.ParallelWriteMM(fname, true, dibella::CkOutputMMHandler());


// fname = name + "-" + std::to_string(nnz) + ".mm";
// nT.ParallelWriteMM(fname, true, dibella::CkOutputMMHandler());

// FullyDistVec<int64_t, dibella::CommonKmers> ReduceV(T.getcommgrid());
// FullyDistVec<int64_t, dibella::CommonKmers> ContigV(T.getcommgrid());

// ContigV =  T.Reduce(Row, ReduceMSR_t(), NullValue);
// ReduceV = nT.Reduce(Row, ReduceMSR_t(), NullValue);

// useExtendedBinOp doesn't seem to be used anywhere, only passed as argument?
// ContigV.EWiseApply(ReduceV, GreaterSR_t(), IsNotEndContigSR_t(), false)