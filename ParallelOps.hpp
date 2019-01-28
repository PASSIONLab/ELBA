//
// Created by Saliya Ekanayake on 12/17/18.
//

#include <chrono>
#include <mpi.h>
#include <memory>

#ifndef LBL_DAL_PARALLEL_OPS_H
#define LBL_DAL_PARALLEL_OPS_H

class ParallelOps {
public:
  static const std::shared_ptr<ParallelOps> init(int *argc, char ***argv);

  void teardown_parallelism();
  ~ParallelOps();

  int world_proc_rank;
  int world_procs_count;

  MPI_Comm MPI_COMM_INSTANCE;

private:
  static std::shared_ptr<ParallelOps> instance;
  ParallelOps(int world_proc_rank, int world_procs_count);
};

#endif //LBL_DAL_PARALLEL_OPS_H

