// Created by Giulia Guidi on 04/02/21.

#include <cmath>
#include <map>
#include <fstream>

#include "TraceUtils.hpp"
#include "kmer/CommonKmers.hpp"
#include "Utils.hpp"
#include "CC.h"

// #define MATRIXPOWER
#define MAXPATHLEN 5000

/*! Namespace declarations */
using namespace combblas;
// typedef ContigSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ContigSRing_t;


// std::vector<std::string> 
// CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& S, std::string& myoutput, TraceUtils tu, 
//    PSpMat<dibella::CommonKmers>::DCCols* spSeq, std::shared_ptr<DistributedFastaData> dfd, int64_t nreads)
// {    

   // float balance = S.LoadImbalance();
   // int64_t nnz   = S.getnnz();

   // std::ostringstream outs;
   // outs.str("");
   // outs.clear();
   // outs << "CreateContig::LoadBalance: " << balance << endl;
   // outs << "CreateContig::nonzeros: "    << nnz     << endl;
   // SpParHelper::Print(outs.str());

   // int64_t nCC = 0;
   // FullyDistVec<int64_t, int64_t> myLabelCC = CC(S, nCC);

   //	std::string myccoutput = myoutput + ".cc";
   //	myLabelCC.ParallelWrite(myccoutput, 1);

   // uint64_t nContig = myLabelCC.Reduce(maximum<int64_t>(), (int64_t)0);
   // nContig++; // because of zero based indexing for cluster

   // std::stringstream ncc;
   // ncc << "nContig: " << nContig << endl;
   // SpParHelper::Print(ncc.str());

   // First4Clust(myLabelCC);
   // HistCC(myLabelCC, nCC);

   // PrintCC(myLabelCC, nCC);

#ifdef MATRIXPOWER

    // PSpMat<dibella::CommonKmers>::MPI_DCCols T = S; // (T) traversal matrix
    // PSpMat<dibella::CommonKmers>::MPI_DCCols ContigM(S.getcommgrid()); // Contig matrix
    // dibella::CommonKmers defaultBVal; 

    // Read vector is gonna multiply the matrix and create contig from there
    // FullyDistVec<int64_t, std::array<char, MAXCONTIGLEN>> ReadVector(S.getcommgrid());
    // char* nt; // NT initial value

    // CustomVectorEntry nt;
    /*
     * A vector of read ids that is my path [0, 1, 3, 67] (paths)
     * A vector of offset [10, 40, 50] (offsets) and offsets.size = paths.size -1 
     * Offset 10 tells me that i have to cut the last 10 bases of read1 and concatenate them to read0, then 40 from read3 and concatenate them to read1 etc. 
    */
    // FullyDistVec<int64_t, char*> ReadVector(S.getcommgrid(), nreads, nt);
    //	auto dcsc = spSeq->GetDCSC();

    // GGGG: fill the read vector with sequences
    // IT * ir; //!< row indices, size nz
    //	for (uint64_t i = 0; i < dcsc->nz; ++i)
    //	{
    //		int64_t lrid = dcsc->ir[i]; // local row idx
    //    seqan::Dna5String rseq = *(dfd->row_seq(lrid));

    //    std::array<char, MAXCONTIGLEN> crseq;
    //    std::string cpprseq;
    //    std::copy(begin(rseq), end(rseq), begin(cpprseq)); // @GGGG: this doesnt work

    //    char* crseq = new char[MAXCONTIGLEN];
    //    crseq = &cpprseq[0]; // C++14

    //    std::cout << rseq << std::endl;
    
    //   ReadVector.SetElement(lrid, crseq); // SetElement only work locally (owner)
    //	}

    // ReadVector.DebugPrint();
    // FullyDistVec<int64_t, std::array<char, MAXSEQLEN>> ContigVector(S.getcommgrid());
    // FullyDistVec<int64_t, char*> ContigVector(S.getcommgrid());

    // do
    // { 
    //     // ContigSR concatenates entries
    //     ReadVector = ContigVector;
    //     ContigVector = SpMV<ContigSR>(S, ReadVector);
                
    // } while (ReadVector != ContigVector); // Once the two vec are identical we're done

    // // GGGG: we know how long are the contig(s) from the CC so we could just extract those or can I use a FullySpDist vector and get only on contig?

    // ContigM.ParallelWriteMM("contig.miracle.mm", true, dibella::CkOutputMMHandler());

#endif

    // std::vector<std::string> myContigSet;
    // @GGGG-TODO: Create contig sequence from connected component matrix
    // {
    //      ...
    // }

    // return myContigSet;
// }

// FullyDistVec<int64_t, dibella::CommonKmers> ReduceV(T.getcommgrid());
// FullyDistVec<int64_t, dibella::CommonKmers> ContigV(T.getcommgrid());

// ContigV =  T.Reduce(Row, ReduceMSR_t(), NullValue);
// ReduceV = nT.Reduce(Row, ReduceMSR_t(), NullValue);

// useExtendedBinOp doesn't seem to be used anywhere, only passed as argument?
// ContigV.EWiseApply(ReduceV, GreaterSR_t(), IsNotEndContigSR_t(), false)
