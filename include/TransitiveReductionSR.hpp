// Created by Giulia Guidi on 09/04/20.

#ifndef __DIBELLA_TR_HPP__
#define __DIBELLA_TR_HPP__

#include "kmer/CommonKmers.hpp"
#include "TraceUtils.hpp"
#include "Utils.hpp"
#include "SR.hpp"

#include <sys/time.h> 
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <cmath>

#include <map>
#include <fstream>

#include <string>
#include <sstream>

/*! Namespace declarations */
using namespace combblas;
using namespace std;

/* TR main function */
void TransitiveReduction(PSpMat<dibella::CommonKmers>::MPI_DCCols& B, TraceUtils tu)
{
    PSpMat<dibella::CommonKmers>::MPI_DCCols BT = B;
    BT.Transpose();
    BT.Apply(OverhangTSRing<dibella::CommonKmers, dibella::CommonKmers>()); 

    if(!(BT == B))
    {
        B += BT;
    }

#ifdef DIBELLA_DEBUG
    tu.print_str("Matrix B += BT: ");
    B.PrintInfo();
    B.ParallelWriteMM("matrixBT.mm", true, dibella::CkOutputMMHandler()); 
#endif

    uint nnz, prev;
    double timeA2 = 0, timeC = 0, timeI = 0, timeA = 0;

    /* Gonna iterate on B until there are no more transitive edges to remove */
    do
    {
        prev = B.getnnz();

        /* Find two-hops neighbors
         * C = B^2
         */
        double start = MPI_Wtime();
        PSpMat<dibella::CommonKmers>::MPI_DCCols F = B;
        PSpMat<dibella::CommonKmers>::MPI_DCCols C = Mult_AnXBn_DoubleBuff<MinPlusSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, F);
        timeA2 += MPI_Wtime() - start;

        C.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);

    #ifdef DIBELLA_DEBUG
        C.ParallelWriteMM("matrixB2.mm", true, dibella::CkOutputMMHandler()); 
        tu.print_str("Matrix C = B^2: ");
        C.PrintInfo();
    #endif
    
        start = MPI_Wtime();
        FullyDistVec<int64_t, dibella::CommonKmers> vA(B.getcommgrid());

        dibella::CommonKmers id; 
        vA = B.Reduce(Row, ReduceMSR_t(), id);
        vA.Apply(PlusFBiSRing<dibella::CommonKmers, dibella::CommonKmers>());

        F.DimApply(Row, vA, Bind2ndSR_t());
        
        timeC += MPI_Wtime() - start;

    #ifdef DIBELLA_DEBUG
        F.ParallelWriteMM("matrixF.mm", true, dibella::CkOutputMMHandler()); 
        tu.print_str("Matrix F = B + FUZZ: ");
        F.PrintInfo();
    #endif

        /* Find transitive edges that can be removed
        * I = F >= C 
        */
        start = MPI_Wtime();
        bool isLogicalNot = false;
        PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(F, C, GreaterBinaryOp<dibella::CommonKmers, dibella::CommonKmers>(), isLogicalNot, id);

        I.Prune(ZeroUnaryOp<bool>(), true);
        
        timeI += MPI_Wtime() - start;
    #ifdef DIBELLA_DEBUG
        I.ParallelWriteMM("matrixI.mm", true, dibella::CkOutputMMHandlerBool());
        tu.print_str("Matrix I = F >= B: ");
        I.PrintInfo();
    #endif

        /* Remove transitive edges
        * B = B .* not(I)
        */ 
        start = MPI_Wtime();
        isLogicalNot = true;
        B = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, I, EWiseMulOp<dibella::CommonKmers, bool>(), isLogicalNot, true);

        /* Prune zero-valued overhang */
        B.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);
        timeA += MPI_Wtime() - start;

    #ifdef DIBELLA_DEBUG
        tu.print_str("Matrix B = B .* not(I): ");
        B.PrintInfo();
    #endif
        nnz = B.getnnz();     
        
    } while (nnz != prev);

    tu.print_str("Matrix B, i.e AAt after transitive reduction: ");
    B.PrintInfo();

 #ifdef DIBELLA_DEBUG
    B.ParallelWriteMM("matrixS.mm", true, dibella::CkOutputMMHandler()); 
    double maxtimeA2, maxtimeC, maxtimeI, maxtimeA;
    
    MPI_Reduce(&timeA2, &maxtimeA2, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeC,  &maxtimeC,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeI,  &maxtimeI,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&timeA,  &maxtimeA,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(myrank == 0)
    {
      std::cout << "TransitiveReduction:TimeA2 = " << maxtimeA2 << std::endl;
      std::cout << "TransitiveReduction:TimeC  = " <<  maxtimeC << std::endl;
      std::cout << "TransitiveReduction:TimeI  = " <<  maxtimeI << std::endl;
      std::cout << "TransitiveReduction:TimeA  = " <<  maxtimeA << std::endl;
    }
 #endif
}

#endif // __DIBELLA_TR_HPP__