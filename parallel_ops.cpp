//
// Created by Saliya Ekanayake on 12/17/18.
//
#include "parallel_ops.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <chrono>

typedef std::chrono::duration<double, std::milli> ms_t;

parallel_ops * parallel_ops::initialize(int *argc, char ***argv) {
  int rank, count;
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &count);
  return new parallel_ops(rank, count);
}

void parallel_ops::teardown_parallelism() {
  MPI_Finalize();
}

parallel_ops::parallel_ops(int world_proc_rank, int world_procs_count) :
    world_proc_rank(world_proc_rank),
    world_procs_count(world_procs_count) {
}

parallel_ops::~parallel_ops() {

}
