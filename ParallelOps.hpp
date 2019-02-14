//
// Created by Saliya Ekanayake on 12/17/18.
//
#ifndef LBL_DAL_PARALLEL_OPS_H
#define LBL_DAL_PARALLEL_OPS_H

#include <chrono>
#include <mpi.h>
#include <memory>
#include <CombBLAS/CombBLAS.h>
#include <CombBLAS/CommGrid.h>

/*! Namespace declarations */
using namespace combblas;

class ParallelOps {
public:
  static const std::shared_ptr<ParallelOps> init(int *argc, char ***argv);

  void teardown_parallelism();
  ~ParallelOps();

  int world_proc_rank;
  int world_procs_count;
  std::shared_ptr<CommGrid> grid;

private:
  static std::shared_ptr<ParallelOps> instance;

  ParallelOps(int world_proc_rank, int world_procs_count);
};

#endif //LBL_DAL_PARALLEL_OPS_H

