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
  if(!flag)
  {
    MPI_Finalize();
  }
}

void ParallelOps::write_file_in_parallel(
        const char* file, const std::string &local_data){

    /*! The following code is adopted from CombBLAS at
     * https://people.eecs.berkeley.edu/~aydin/CombBLAS/html/_fully_dist_sp_vec_8cpp_source.html#l01310.
     */
    auto *bytes = new int64_t[world_procs_count];
    bytes[world_proc_rank] = local_data.size();
    MPI_Allgather(MPI_IN_PLACE, 1, MPIType<int64_t>(), bytes, 1,
                  MPIType<int64_t>(), MPI_COMM_WORLD);
    /*! Note. std::accumulate has an open range on the right meaning
     * that, for example, the 'bytesuntil' is the sum of bytes upto my
     * rank and not including my bytes.
     */
    int64_t bytesuntil = std::accumulate(bytes, bytes + world_proc_rank,
                                         static_cast<int64_t>(0));
    int64_t bytestotal = std::accumulate(bytes, bytes + world_procs_count,
                                         static_cast<int64_t>(0));

    if (world_proc_rank == 0) {
        // only leader rights the original file with no content
        std::ofstream ofs(file, std::ios::out);
#ifndef NDEBUG
        std::cout << "Creating file" << file << " with "
        << bytestotal << " bytes" << std::endl;
#endif
        ofs.seekp(bytestotal - 1);
        // this will likely create a sparse file so the actual disks won't spin yet
        ofs.write("", 1);
        ofs.close();
    }
    MPI_Barrier(MPI_COMM_WORLD);

    struct stat st;     // get file row_size
    if (stat(file, &st) == -1) {
        MPI_Abort(MPI_COMM_WORLD, NOFILE);
    }
    if (world_proc_rank ==
        world_procs_count - 1)  // let some other processor do the testing
    {
#ifndef NDEBUG
        std::cout << "File " << file << " is actually " << st.st_size
                  << " bytes seen from process " << world_proc_rank
                  << std::endl;
#endif
    }

    // Then everyone fills it
    FILE *ffinal;
    if ((ffinal = fopen(file, "r+")) == NULL) {
        printf("ERROR: Output file %s failed to open at process %d\n",
               file, world_proc_rank);
        MPI_Abort(MPI_COMM_WORLD, NOFILE);
    }

    fseek(ffinal, bytesuntil, SEEK_SET);
    fwrite(local_data.c_str(), 1, bytes[world_proc_rank], ffinal);
    fflush(ffinal);
    fclose(ffinal);
    delete[] bytes;
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