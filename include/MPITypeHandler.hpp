#ifndef MPI_TYPE_MAKER_H_
#define MPI_TYPE_MAKER_H_

#include "common.h"
#include <mpi.h>
#include <iostream>

struct MPITypeHandler
{
    MPI_Datatype *datatype;

    MPITypeHandler(MPI_Datatype *dtype) : datatype(dtype) { MPI_Type_commit(datatype); }
    MPI_Datatype getdtype() const { return *datatype; }

    ~MPITypeHandler() { MPI_Type_free(datatype); }
};


#endif
