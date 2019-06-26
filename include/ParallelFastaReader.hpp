// Created by Saliya Ekanayake on 1/6/19.

#ifndef LBL_DAL_PARALLEL_FASTA_READER_HPP
#define LBL_DAL_PARALLEL_FASTA_READER_HPP


#include "ParallelOps.hpp"
#include "FastaData.hpp"
#include "TraceUtils.hpp"
#include "DistributedFastaData.hpp"


class ParallelFastaReader {
public:
  static void read_fasta(const char *file, uint64_t overlap, int rank,
                  int world_size, char *&buff,
                  uint64_t &l_start, uint64_t &l_end);
};


#endif //LBL_DAL_PARALLEL_FASTA_READER_HPP
