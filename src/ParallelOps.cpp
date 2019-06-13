// Created by Saliya Ekanayake on 12/17/18.

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include "../include/ParallelOps.hpp"

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
  int flag;
  MPI_Initialized(&flag);
  if(!flag) {
    MPI_Finalize();
  }
}

ParallelOps::ParallelOps(int world_proc_rank, int world_procs_count) :
    world_proc_rank(world_proc_rank),
    world_procs_count(world_procs_count) {

  auto grid_size = static_cast<int>(sqrt(world_procs_count));
  grid = std::make_shared<CommGrid>(MPI_COMM_WORLD, grid_size, grid_size);
}

ParallelOps::~ParallelOps() {
  int flag;
  MPI_Initialized(&flag);
  if (!flag) {
    MPI_Finalize();
  }
}