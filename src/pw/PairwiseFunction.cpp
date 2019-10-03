// Created by Saliya Ekanayake on 2019-07-10.

#include "../../include/pw/PairwiseFunction.hpp"


PairwiseFunction::PairwiseFunction()= default;

PairwiseFunction::~PairwiseFunction() = default;

void PairwiseFunction::add_time(std::string type, double duration) {
  if (types.find(type) != types.end()){
    size_t idx = types[type];
    ++counts[idx];
    times[idx] += duration;
  } else {
    types[type] = times.size();
    counts.push_back(1);
    times.push_back(duration);
  }
}

void PairwiseFunction::print_avg_times(std::shared_ptr<ParallelOps> parops) {
  MPI_Allreduce(MPI_IN_PLACE, &(counts[0]), counts.size(),
      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &(times[0]), times.size(),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(parops->world_proc_rank == 0) {
    for (auto &type : types) {
      std::cout << "  PWF:" << type.first << " avg per proc:"
                << (times[type.second] / parops->world_procs_count) << " ms" << std::endl;
    }
  }
}