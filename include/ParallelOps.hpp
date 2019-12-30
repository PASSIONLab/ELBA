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
  void write_file_in_parallel(const char* file, const std::string &local_data);
  std::string bcast(std::string str){
    if(world_proc_rank == 0) {
      MPI_Bcast(&str[0], str.length(), MPI_CHAR, 0, MPI_COMM_WORLD);
      return str;
    } else {
      char* buf = new char[str.length()];
      MPI_Bcast(buf, str.length(), MPI_CHAR, 0, MPI_COMM_WORLD);
      std::string s(buf);
      delete [] buf;
      return s;
    }
  }
  ~ParallelOps();

  int world_proc_rank;
  int world_procs_count;
  std::shared_ptr<CommGrid> grid;

private:
  static std::shared_ptr<ParallelOps> instance;

  ParallelOps(int world_proc_rank, int world_procs_count);
};

#endif //LBL_DAL_PARALLEL_OPS_H

