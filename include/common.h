#ifndef COMMON_H_
#define COMMON_H_

#include <vector>
#include <array>
#include <string>
#include <tuple>
#include <set>
#include <unordered_map>
#include <memory>
#include <cstdint>
#include <limits>
#include <mpi.h>
#include "CombBLAS/CombBLAS.h"

using namespace combblas;
using Grid = std::shared_ptr<CommGrid>;

#if MPI_VERSION == 3
#define MPI_HAS_LARGE_COUNTS 0

using MPI_Count_type = int;
using MPI_Displ_type = int;

#define MPI_COUNT_TYPE MPI_INT
#define MPI_FUNC_SELECT(funcname) funcname

#elif MPI_VERSION == 4
#define MPI_HAS_LARGE_COUNTS 1

using MPI_Count_type = MPI_Count;
using MPI_Displ_type = MPI_Aint;

#define MPI_COUNT_TYPE MPI_COUNT
#define MPI_FUNC_SELECT(funcname) funcname##_c

#else
#error "MPI version should either be 3 or 4."
#endif

#define MPI_ALLTOALL         MPI_FUNC_SELECT(MPI_Alltoall)
#define MPI_ALLTOALLV        MPI_FUNC_SELECT(MPI_Alltoallv)
#define MPI_ALLREDUCE        MPI_FUNC_SELECT(MPI_Allreduce)
#define MPI_SCATTER          MPI_FUNC_SELECT(MPI_Scatter)
#define MPI_SCATTERV         MPI_FUNC_SELECT(MPI_Scatterv)
#define MPI_GATHER           MPI_FUNC_SELECT(MPI_Gather)
#define MPI_GATHERV          MPI_FUNC_SELECT(MPI_Gatherv)
#define MPI_ALLGATHER        MPI_FUNC_SELECT(MPI_Allgather)
#define MPI_ALLGATHERV       MPI_FUNC_SELECT(MPI_Allgatherv)
#define MPI_BCAST            MPI_FUNC_SELECT(MPI_Bcast)
#define MPI_REDUCE           MPI_FUNC_SELECT(MPI_Reduce)
#define MPI_ISEND            MPI_FUNC_SELECT(MPI_Isend)
#define MPI_IRECV            MPI_FUNC_SELECT(MPI_Irecv)
#define MPI_FILE_READ_AT_ALL MPI_FUNC_SELECT(MPI_File_read_at_all)

#ifndef MPI_SIZE_T
static_assert(std::numeric_limits<size_t>::max() == std::numeric_limits<unsigned long>::max());
#define MPI_SIZE_T MPI_UNSIGNED_LONG
#endif

template <class NT>
struct CT
{
    typedef SpDCCols<uint64_t, NT> PSpDCCols;
    typedef SpParMat<uint64_t, NT, PSpDCCols> PSpParMat;
    typedef FullyDistVec<uint64_t, NT> PDistVec;
};

#ifndef USE_BLOOM
#define USE_BLOOM 1
#endif

#endif
