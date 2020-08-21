// Created by Saliya Ekanayake on 2019-02-13.

#include "../include/TraceUtils.hpp"

void TraceUtils::print_msg(const std::string &title, const std::string &msg,
                           const std::shared_ptr<ParallelOps> &parops) {
  int flag;
  if (parops->world_proc_rank > 0)
    MPI_Recv(&flag, 1, MPI_INT, parops->world_proc_rank - 1,
             99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::cout << "\nRank: " << parops->world_proc_rank << " - " << title
            << std::endl;
  std::cout << msg << std::endl;
  if (parops->world_proc_rank < parops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, parops->world_proc_rank + 1, 99,
             MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void
TraceUtils::print_msg_on_rank(const std::string &title, const std::string &msg,
                              const std::shared_ptr<ParallelOps> &parops,
                              int rank) {
  if (parops->world_proc_rank == rank) {
    std::cout << "\nRank: " << parops->world_proc_rank << " - " << title
              << std::endl;
    std::cout << msg << std::endl;
  }

}

TraceUtils::TraceUtils(bool is_print_rank) : is_print_rank(is_print_rank) {

}

void TraceUtils::print_str(std::string str) {
  if(is_print_rank){
    std::cout<<str;
    std::cout.flush();
  }
}


