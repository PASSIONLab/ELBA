// Created by Saliya Ekanayake on 2019-07-05.

#ifndef LBL_PISA_PAIRWISEFUNCTION_HPP
#define LBL_PISA_PAIRWISEFUNCTION_HPP

#include <unordered_map>
#include <string>
#include <seqan/score.h>
#include "../Kmer.hpp"
#include "../AlignmentInfo.hpp"
#include "../kmer/CommonKmers.hpp"

class PairwiseFunction {
public:
  PairwiseFunction();
  virtual ~PairwiseFunction();

  virtual void apply(uint64_t l_col_idx, uint64_t g_col_idx,
      uint64_t l_row_idx, uint64_t g_row_idx,
      seqan::Peptide *seq_h, seqan::Peptide *seq_v,
      pisa::CommonKmers &cks, std::stringstream& ss) = 0;

  void add_time(std::string type, double duration);
  void print_avg_times(std::shared_ptr<ParallelOps> parops);

private:
  std::unordered_map<std::string, size_t> types;
  std::vector<uint64_t> counts;
  std::vector<double> times;

};

#endif //LBL_PISA_PAIRWISEFUNCTION_HPP
