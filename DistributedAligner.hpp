//
// Created by Saliya Ekanayake on 2019-02-19.
//

#ifndef LBL_DAL_ALIGNER_H
#define LBL_DAL_ALIGNER_H

#include <seqan/align.h>
#include "Utils.hpp"
#include "Kmer.hpp"
#include "DistributedFastaData.h"

class DistributedAligner {
public:
  DistributedAligner(ushort seed_length, int xdrop, int gap_open, int gap_ext, const std::shared_ptr<DistributedFastaData> dfd,
                     PSpMat<CommonKmers>::MPI_DCCols mat,
                     const std::shared_ptr<ParallelOps> &parops);

  void align();

private:
  ushort seed_length;
  int xdrop;
  int gap_open;
  int gap_ext;
  PSpMat<CommonKmers>::MPI_DCCols mat;
  std::shared_ptr<DistributedFastaData> dfd;
  std::shared_ptr<ParallelOps> parops;

  std::vector<seqan::AlignmentStats> scores;

};


#endif //LBL_DAL_ALIGNER_H
