// #include <cmath>
// #include "../include/Constants.hpp"
// #include "../include/ParallelOps.hpp"
// #include "../include/ParallelFastaReader.hpp"
// #include "../include/Alphabet.hpp"
// #include "../include/DistributedPairwiseRunner.hpp"
// #include "../include/cxxopts.hpp"
// #include "../include/TransitiveReductionSR.hpp"
// #include "../include/Utils.hpp"
// #include "../include/CC.h"

// #include <map>
// #include <fstream>

// /*! Namespace declarations */
// using namespace combblas;

// typedef BestOverlapVSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> BestOvVSR_t;
// typedef BestOverlapMSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> BestOvMSR_t;

// std::vector<std::string> 
// CreateContig(PSpMat<dibella::CommonKmers>::MPI_DCCols& B, std::string& myoutput, TraceUtils tu)
// {
//     /* (1) Find min value (best/longest overlap) per column entry and store it into a vector */
//     FullyDistVec<int64_t, dibella::CommonKmers> vB(B.getcommgrid());
//     dibella::CommonKmers id;
//     vB = B.Reduce(Column, BestOvVSR_t(), id);

//     /* Best overlap matrix (set entry to zero if != bestOverlapV entry) 
//      * (2) DimApply and remove if entry != vector entry
//      */
//     B.DimApply(Column, vB, BestOvMSR_t());

//     /* Prune zero-valued overhang entries */
//     B.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);
//     B.PrintInfo();

//     /* (3) Connected component for contig extraction (from LACC: https://github.com/PASSIONLab/CombBLAS/blob/e6c55bd48a442b8fa95870fb5f18cd6b89cbffe9/Applications/CC.cpp)
//      * 
//      * The matrix should be already symmetric at this point.    
//      * If they need to be identical I can save both and then agree than row > col must read only [0] and col > row only [1]?
//      * 
     
//     BT.Transpose();
//     if(!(BT == B))
//     {
//         SpParHelper::Print("Symmatricizing an unsymmetric input matrix.\n");
//         B += BT;
//     }
//     B.PrintInfo();

//     *
//     */

//     float balance = B.LoadImbalance();
//     int64_t nnz   = B.getnnz();

//     std::ostringstream outs;
//     outs.str("");
//     outs.clear();
//     outs << "CreateContig::LoadBalance: " << balance << endl;
//     outs << "CreateContig::nonzeros: "    << nnz     << endl;
//     SpParHelper::Print(outs.str());

//     uint64_t nCC = 0;
//     FullyDistVec<uint64_t, uint64_t> myLabelCC = CC(B, nCC);

// 	std::string myccoutput = myoutput + ".cc";
// 	myLabelCC.ParallelWrite(myccoutput, 1);

//     uint64_t nContig = myLabelCC.Reduce(maximum<uint64_t>(), (uint64_t)0);
//     nContig++; // because of zero based indexing for cluster

//     std::stringstream ncc;
//     ncc << "nContig: " << nContig << endl;
//     SpParHelper::Print(ncc.str());

//     tu.print_str("BestOverlap matrix: ");
//     B.PrintInfo();

//     std::vector<std::string> myContigSet;
//     // @GGGG-TODO: Create contig sequence from connected component matrix
//     // {
//     //      ...
//     // }

//     return myContigSet;
// }