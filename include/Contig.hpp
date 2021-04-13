// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "SR.hpp"
#include "CC.h"

#define MATRIXPOWER

/*! Namespace declarations */
using namespace combblas;
typedef ContigSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ContigSRing_t;

std::vector<std::string> 
CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu)
{    

#ifdef CC
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

    First4Clust(myLabelCC);
    HistCC(myLabelCC, nCC);
    PrintCC(myLabelCC, nCC);
#endif

#ifdef MATRIXPOWER

    uint nnz;

    PSpMat<dibella::CommonKmers>::MPI_DCCols T = S; // (T) traversal matrix
    PSpMat<dibella::CommonKmers>::MPI_DCCols ContigM(S.getcommgrid()); // Contig matrix

    dibella::CommonKmers defaultBVal; 

    do
    { 
        // if T(i,j) is end of contig do nothing (return id)
        PSpMat<dibella::CommonKmers>::MPI_DCCols P = // (P) power matrix
            Mult_AnXBn_DoubleBuff<ContigSRing_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(S, T);

        // if S(i,j).cend == true, then P(i,j).cend == true
        P = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(P, S, 
            CheckEndContigSR<dibella::CommonKmers, dibella::CommonKmers>(), false, defaultBVal);

        // if P(i,j).cend == true, copy it over since it's a complete contig
        ContigM = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(ContigM, P, 
            CopyOverContigSR<dibella::CommonKmers, dibella::CommonKmers>(), false, defaultBVal);
            
        P.Prune(ZeroOverhangContigSR<dibella::CommonKmers>(), true);

        // update matrix (if cend = true, treat nonzero as zero and move on)
        T = P;
        nnz = P.getnnz();
                
    } while (nnz != 0); // i might need a new termination condition

    ContigM.ParallelWriteMM("contig.miracle.mm", true, dibella::CkOutputMMHandler());

#endif

    std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    return myContigSet;
}

// FullyDistVec<int64_t, dibella::CommonKmers> ReduceV(T.getcommgrid());
// FullyDistVec<int64_t, dibella::CommonKmers> ContigV(T.getcommgrid());

// ContigV =  T.Reduce(Row, ReduceMSR_t(), NullValue);
// ReduceV = nT.Reduce(Row, ReduceMSR_t(), NullValue);

// useExtendedBinOp doesn't seem to be used anywhere, only passed as argument?
// ContigV.EWiseApply(ReduceV, GreaterSR_t(), IsNotEndContigSR_t(), false)