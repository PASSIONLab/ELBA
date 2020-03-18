// Created by Saliya Ekanayake on 2019-07-10.

#include <algorithm>

#include "../../include/pw/PairwiseFunction.hpp"


// PairwiseFunction::PairwiseFunction()= default;
PairwiseFunction::PairwiseFunction() :
	nalignments(0)
{
}

PairwiseFunction::~PairwiseFunction() = default;

void PairwiseFunction::add_time(std::string type, double duration) {
  // may be called threaded
  int tid = 0;
  #ifdef THREADED
  tid = omp_get_thread_num();
  #endif
  
  if (types[tid].find(type) != types[tid].end()){
    size_t idx = types[tid][type];
    ++counts[tid][idx];
    times[tid][idx] += duration;
  } else {
	types[tid][type] = times[tid].size();
    counts[tid].push_back(1);
    times[tid].push_back(duration);
  }
}

void PairwiseFunction::print_avg_times(std::shared_ptr<ParallelOps> parops) {
  // counts and times can be empty if there were non non-zeros in that
  // process grid cell. Therefore, we need to fix the collective call as follows.
  int counts_size = 0;
  for (auto el : counts)
	counts_size = std::max(counts_size, static_cast<int>(el.size()));

  MPI_Allreduce(MPI_IN_PLACE, &counts_size, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // std::cout << "counts_size after allreduce " << counts_size << std::endl;
  // for (int tid = 0; tid < MAX_THD; ++tid)
  // {
  // 	  if (types[tid].size() != 0)
  // 	  {
  // 		  std::cout << "tid " << tid << " stats:\n";
  // 		  for (auto &type : types[tid])
  // 		  {
  // 			  std::cout << "  PWF: "  << type.first << " - time: "
  // 						<< times[tid][type.second] << " - count: "
  // 						<< counts[tid][type.second] << "\n";
  // 		  }
  // 	  }
  // }

  // max over all threads
  std::vector<uint64_t> counts_tmp(counts_size, 0);  
  for (auto &el : counts)
    for (int i = 0; i < el.size(); ++i)
	  counts_tmp[i] = std::max(counts_tmp[i], el[i]);
  
  std::vector<double> times_tmp(counts_size, 0);
  for (auto &el : times)
    for (int i = 0; i < el.size(); ++i)
	  times_tmp[i] = std::max(times_tmp[i], el[i]);
  
  std::unordered_map<std::string, size_t> *types_tmp = NULL;
  for (auto &el : types)
  {
	  if (el.size() != 0)
	  {
		  types_tmp = &el;
		  break;
	  }
  }

  MPI_Allreduce(MPI_IN_PLACE, &(counts_tmp[0]), counts_tmp.size(),
      MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE, &(times_tmp[0]), times_tmp.size(),
                MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (parops->world_proc_rank == 0 && types_tmp != NULL)
  {
    for (auto &type : *types_tmp)
	{
      std::cout << "  PWF:" << type.first << " avg per proc:"
                << (times_tmp[type.second] / parops->world_procs_count)
				<< " ms" << std::endl;
    }
  }
}
