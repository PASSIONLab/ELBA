//
// Created by Saliya Ekanayake on 1/6/19.
//

#ifndef LBL_DAL_PARALLEL_FASTA_READER_HPP
#define LBL_DAL_PARALLEL_FASTA_READER_HPP


#include "ParallelOps.h"

class ParallelFastaReader {
public:
  void readFasta(const char *file, int seq_count, int overlap, int rank, int world_size);
};


#endif //LBL_DAL_PARALLEL_FASTA_READER_HPP
