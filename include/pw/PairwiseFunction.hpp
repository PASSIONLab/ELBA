// Created by Saliya Ekanayake on 2019-07-05.

#ifndef DIBELLA_PAIRWISEFUNCTION_HPP
#define DIBELLA_PAIRWISEFUNCTION_HPP

#include <unordered_map>
#include <string>
#include <seqan/score.h>
#include <seqan/align_parallel.h>
#include "../AlignmentInfo.hpp"
#include "../kmer/CommonKmers.hpp"
#include "../ParallelOps.hpp"
#include "../DistributedFastaData.hpp"
#include "../Utils.hpp"

class PairwiseFunction {
public:

  static const int MAX_THD = 128;
	
  PairwiseFunction();
  virtual ~PairwiseFunction();

  virtual void apply(uint64_t l_col_idx, uint64_t g_col_idx,
      uint64_t l_row_idx, uint64_t g_row_idx,
      seqan::Dna5String *seq_h, seqan::Dna5String *seq_v,
      ushort k,
      elba::CommonKmers &cks, std::stringstream& ss) = 0;

  virtual
  void
  apply_batch (
         seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsh,
			   seqan::StringSet<seqan::Gaps<seqan::Dna5String>> &seqsv,
			   uint64_t *lids,
			   uint64_t col_offset,
			   uint64_t row_offset,
			   PSpMat<elba::CommonKmers>::ref_tuples *mattuples,
         std::ofstream &lfs,
         const bool noAlign,
         ushort k,
         uint64_t nreads,
         std::vector<int64_t>& ContainedSeqPerProc,
         float ratioScoreOverlap = 0.99,     // GGGG: Precomputed for error rate = 15% and default scoring matrix (1,-1,-1) (0.445 for CLR, 0.99 for CCS)
         //float ratioScoreOverlap = 0.445,     // GGGG: Precomputed for error rate = 15% and default scoring matrix (1,-1,-1) (0.445 for CLR, 0.99 for CCS)
	 int debugThr = 50) = 0;            // GGGG: Fixed threshold, this is convenient only for debugging


  void add_time(std::string type, double duration);
  void print_avg_times(std::shared_ptr<ParallelOps> parops, std::ofstream &lfs);

  uint64_t nalignments;
  
private:
  std::unordered_map<std::string, size_t>	types[MAX_THD];
  std::vector<uint64_t>	counts[MAX_THD];
  std::vector<double>		times[MAX_THD];
};

#endif //DIBELLA_PAIRWISEFUNCTION_HPP
