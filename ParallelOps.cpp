//
// Created by Saliya Ekanayake on 12/17/18.
//
#include "ParallelOps.hpp"
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <chrono>

typedef std::chrono::duration<double, std::milli> ms_t;

std::shared_ptr<ParallelOps> ParallelOps::instance = nullptr;

const std::shared_ptr<ParallelOps> ParallelOps::init(int *argc, char ***argv) {
  if (instance != nullptr)
    return instance;

  int rank, count;
  MPI_Init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &count);
  instance = std::shared_ptr<ParallelOps>(new ParallelOps(rank, count));
  return instance;
}

void ParallelOps::teardown_parallelism() {
  MPI_Finalize();
}

ParallelOps::ParallelOps(int world_proc_rank, int world_procs_count) :
    world_proc_rank(world_proc_rank),
    world_procs_count(world_procs_count) {
}

ParallelOps::~ParallelOps() {
  int *flag = nullptr;
  MPI_Initialized(flag);
  if (flag){
    MPI_Finalize();
  }
}
