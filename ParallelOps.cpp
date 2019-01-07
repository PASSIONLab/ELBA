//
// Created by Saliya Ekanayake on 12/17/18.
//
#include "ParallelOps.h"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <chrono>

typedef std::chrono::duration<double, std::milli> ms_t;

ParallelOps * ParallelOps::initialize(int *argc, char ***argv) {
  int rank, count;
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &count);
  return new ParallelOps(rank, count);
}

void ParallelOps::teardown_parallelism() {
  MPI_Finalize();
}

ParallelOps::ParallelOps(int world_proc_rank, int world_procs_count) :
    world_proc_rank(world_proc_rank),
    world_procs_count(world_procs_count) {
}

ParallelOps::~ParallelOps() {

}
