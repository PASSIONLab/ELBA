// Created by Saliya Ekanayake on 2/10/19.

#ifndef LBL_DAL_DEBUGUTILS_HPP
#define LBL_DAL_DEBUGUTILS_HPP

#include "FastaData.hpp"
#include "ParallelOps.hpp"

void _debug_print_fasta_data(std::shared_ptr<FastaData> fd,
                      std::shared_ptr<ParallelOps> p_ops){
  /* Debug print of FastData instance*/
  int flag;
  if (p_ops->world_proc_rank > 0)
    MPI_Recv(&flag, 1, MPI_INT, p_ops->world_proc_rank-1,
             99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::cout<<"\nRank: "<<p_ops->world_proc_rank<<"\n------------------------\n";
  fd->print();
  if (p_ops->world_proc_rank < p_ops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, p_ops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

#endif //LBL_DAL_DEBUGUTILS_HPP
