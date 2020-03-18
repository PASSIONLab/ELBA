// Created by Saliya Ekanayake on 2019-07-05.

#ifndef DISTAL_PAIRWISEFUNCTION_HPP
#define DISTAL_PAIRWISEFUNCTION_HPP

#include <unordered_map>
#include <string>
#include <seqan/score.h>
#include "../AlignmentInfo.hpp"
#include "../kmer/CommonKmers.hpp"
#include "../ParallelOps.hpp"

class PairwiseFunction {
public:

  static const int MAX_THD = 128;
	
  PairwiseFunction();
  virtual ~PairwiseFunction();

  virtual void apply(uint64_t l_col_idx, uint64_t g_col_idx,
      uint64_t l_row_idx, uint64_t g_row_idx,
      seqan::Peptide *seq_h, seqan::Peptide *seq_v,
      distal::CommonKmers &cks, std::stringstream& ss) = 0;

  void add_time(std::string type, double duration);
  void print_avg_times(std::shared_ptr<ParallelOps> parops);

  uint64_t nalignments;
  
private:
  std::unordered_map<std::string, size_t>	types[MAX_THD];
  std::vector<uint64_t>						counts[MAX_THD];
  std::vector<double>						times[MAX_THD];
};

#endif //DISTAL_PAIRWISEFUNCTION_HPP
