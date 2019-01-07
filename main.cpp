#include <iostream>
//#include <boost/program_options/options_description.hpp>
//#include <boost/program_options/option.hpp>
//#include <boost/program_options/variables_map.hpp>
//#include <boost/program_options/parsers.hpp>
#include <cxxopts.hpp>
#include "constants.h"
#include "ParallelOps.h"
#include "ParallelFastaReader.hpp"

/* Namespace declarations */
//namespace po = boost::program_options;

/* Function signatures */
int parse_args(int argc, char **argv);

/* Global variables */
ParallelOps *p_ops;
std::string input_file;
int input_overlap;
int input_seq_count;

int main(int argc, char **argv) {
  p_ops = ParallelOps::initialize(&argc, &argv);
  int ret = parse_args(argc, argv);
  if (ret < 0){
    p_ops->teardown_parallelism();
    return ret;
  }

  ParallelFastaReader pfr;
  pfr.readFasta(input_file.c_str(), input_seq_count, input_overlap, p_ops->world_proc_rank, p_ops->world_procs_count);

//  std::cout << "Hello, World!" << std::endl;
  return 0;
}

int parse_args(int argc, char **argv) {
  cxxopts::Options options("Distributed Aligner", "A distributed protein aligner");

  options.add_options()
      (CMD_OPTION_INPUT, CMD_OPTION_DESCRIPTION_INPUT, cxxopts::value<std::string>())
      (CMD_OPTION_INPUT_SEQ_COUNT, CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT, cxxopts::value<int>())
      (CMD_OPTION_INPUT_OVERLAP, CMD_OPTION_DESCRIPTION_INPUT_OVERLAP, cxxopts::value<int>())
      ;

  auto result = options.parse(argc, argv);

  bool is_world_rank0 = p_ops->world_proc_rank == 0;
  if (result.count(CMD_OPTION_INPUT)){
    input_file = result[CMD_OPTION_INPUT].as<std::string>();
  }else {
    if (is_world_rank0)
      std::cout<<"ERROR: Input file not specified"<<std::endl;
    return -1;
  }

  if (result.count(CMD_OPTION_INPUT_SEQ_COUNT)){
    input_seq_count = result[CMD_OPTION_INPUT_SEQ_COUNT].as<int>();
  }else {
    if (is_world_rank0)
      std::cout<<"ERROR: Input sequence count not specified"<<std::endl;
    return -1;
  }

  if (result.count(CMD_OPTION_INPUT_OVERLAP)){
    input_overlap = result[CMD_OPTION_INPUT_OVERLAP].as<int>();
  }else {
    input_overlap = 10000;
  }

  return 0;
}
