// Created by Giulia Guidi on 09/04/20.

#include "CombBLAS/CombBLAS.h"
#include "../include/TransitiveReductionSR.hpp"
#include "../include/Defines.hpp"
#include <mpi.h>
#include <sys/time.h> 
#include <iostream>
#include <functional>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

int cblas_splits;

double cblas_alltoalltime;
double cblas_allgathertime;
double cblas_mergeconttime;
double cblas_transvectime;
double cblas_localspmvtime;

/** Given a biridrected graph, an edge v ?-? x can only be considered transitive given a pair of edges v ?-? w ?-? x if:
 * (1) The two heads adjacent to w have opposite orientation: 
 *      2nd bit != 1st bit in MinPlus semiring B = A^2 such as 01 and 01 or 10 and 10;
 * (2) The heads adjacent to v in v ?-? w and v ?-? x have the same orientation, and
 * (3) The heads adjacent to x in v ?-? x and w ?-? x have the same orientation:
 *      1st and 2nd bit in M == 1st and 2nd bit in B during I = M >= B.
*/

// GGGG: TODO add fields to dibella::CommonKmers
PSpMat<dibella::CommonKmers>::MPI_DCCols
PerformTransitiveReduction(PSpMat<dibella::CommonKmers>::MPI_DCCols& A)
{
	int nprocs, myrank;
	int cblas_splits = 1;	
    
	MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);

	if(myrank == 0)
	{
	        cout << "Hey, transitive reduction is starting"   << endl;
	}

        typedef MinPlusBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> MinPlusSR;
        typedef Bind2ndBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> Bind2ndSR;
        typedef ReduceMBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ReduceMSR;
        
        dibella::CommonKmers id;

	// shared_ptr<CommGrid> fullWorld;
	// fullWorld.reset(new CommGrid(MPI_COMM_WORLD, 0, 0));
        // PSpMat<dibella::CommonKmers>::MPI_DCCols A(fullWorld);

#ifdef DEBUG
        A.PrintInfo();
#endif

        /* Randomly permute for load balance */
        // int perm = false; // GGGG: TODO make this input parameter
        // if(perm == 1)
        // {
        //     SpParHelper::Print("Performing random permutation of matrix\n");

        //     FullyDistVec<int64_t, dibella::CommonKmers> prow(A.getcommgrid());
        //     FullyDistVec<int64_t, dibella::CommonKmers> pcol(A.getcommgrid());

        //         // GGGG: from line 931 of FullyDistVect.cpp:
        //         // ABAB: In its current form, unless LengthUntil returns NT
        //         // this won't work reliably for anything other than integers
        //         // GGGG: TODO look into this
        //
        //     prow.iota(A.getnrow(), id);
        //     pcol.iota(A.getncol(), id);

        //     prow.RandPerm();
        //     pcol.RandPerm();
        //     (A)(prow, pcol, true);

        //     SpParHelper::Print("Performed random permutation of matrix\n");
        // }

        // GGGG: B = A^2
        PSpMat<dibella::CommonKmers>::MPI_DCCols M = A;
        PSpMat<dibella::CommonKmers>::MPI_DCCols B = Mult_AnXBn_DoubleBuff<MinPlusSR, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(A, M);
#ifdef DEBUG
        B.PrintInfo();
#endif

//         FullyDistVec<int64_t, dibella::CommonKmers> vA(M.getcommgrid());
        
//         // GGGG: give meaning to dibella::CommonKmers());
//         vA = M.Reduce(Row, ReduceMSR(), id);
//         vA.Apply(PlusFBiSRing<dibella::CommonKmers, dibella::CommonKmers>());
// #ifdef DEBUG
//         vA.DebugPrint();
// #endif

//         M.DimApply(Row, vA, Bind2ndSR());
// #ifdef DEBUG
//         M.PrintInfo();
// #endif

//         // GGGG: I = M >= B 
//         bool isLogicalNot = false;
//         PSpMat<bool>::MPI_DCCols I = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(M, B, GreaterBinaryOp<dibella::CommonKmers, dibella::CommonKmers>(), isLogicalNot, id);

//         // GGGG: prune potential zero-valued nonzeros
//         I.Prune(ZeroUnaryOp<dibella::CommonKmers>(), true);
// #ifdef DEBUG
//         I.PrintInfo();
// #endif

//         // GGGG: A = A .* not(I)
//         isLogicalNot = true;
//         PSpMat<dibella::CommonKmers>::MPI_DCCols C = EWiseMult(A, I, isLogicalNot);
// #ifdef DEBUG
//         C.PrintInfo();
// #endif
// 	return C;
}
