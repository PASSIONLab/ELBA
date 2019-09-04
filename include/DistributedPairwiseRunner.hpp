//
// Created by Saliya Ekanayake on 2019-02-19.
//

#ifndef LBL_DAL_ALIGNER_H
#define LBL_DAL_ALIGNER_H

#include <seqan/align.h>
#include "Utils.hpp"
#include "Kmer.hpp"
#include "DistributedFastaData.hpp"
#include "pw/PairwiseFunction.hpp"
#include "AlignmentInfo.hpp"

//struct AlignmentInfo{
//  seqan::AlignmentStats stats;
//  ushort seq_h_length;
//  ushort seq_v_length;
//  ushort seq_h_seed_length;
//  ushort seq_v_seed_length;
//  uint64_t seq_h_g_idx;
//  uint64_t seq_v_g_idx;
//};

class DistributedPairwiseRunner {
public:
  DistributedPairwiseRunner(std::shared_ptr<DistributedFastaData> dfd,
                     PSpMat<CommonKmers>::MPI_DCCols mat,
                     const std::shared_ptr<ParallelOps> &parops);

//  uint64_t align_seqs();
  void write_overlaps(const char *file);
  void run(PairwiseFunction *pf);

private:
  PSpMat<CommonKmers>::MPI_DCCols mat;
  std::shared_ptr<DistributedFastaData> dfd;
  std::shared_ptr<ParallelOps> parops;
};


#endif //LBL_DAL_ALIGNER_H
