#include <iostream>
//#include <boost/program_options/options_description.hpp>
//#include <boost/program_options/option.hpp>
//#include <boost/program_options/variables_map.hpp>
//#include <boost/program_options/parsers.hpp>
#include <cxxopts.hpp>
#include "constants.h"
#include "parallel_ops.h"

/* Namespace declarations */
//namespace po = boost::program_options;

/* Function signatures */
int parse_args(int argc, char **argv);

/* Global variables */
parallel_ops *p_ops;
std::string input_file;

int main(int argc, char **argv) {
  p_ops = parallel_ops::initialize(&argc, &argv);
  int ret = parse_args(argc, argv);
  if (ret < 0){
    p_ops->teardown_parallelism();
    return ret;
  }


  std::cout << "Hello, World!" << std::endl;
  return 0;
}

int parse_args(int argc, char **argv) {
  cxxopts::Options options("Distributed Aligner", "A distributed protein aligner");

  options.add_options()
      (CMD_OPTION_SHORT_INPUT, CMD_OPTION_DESCRIPTION_INPUT, cxxopts::value<std::string>())
      ;

  auto result = options.parse(argc, argv);

  bool is_world_rank0 = p_ops->world_proc_rank == 0;
  if (result.count(CMD_OPTION_SHORT_INPUT)){
    input_file = result[CMD_OPTION_SHORT_INPUT].as<std::string>();
  }else {
    if (is_world_rank0)
      std::cout<<"ERROR: Input file not specified"<<std::endl;
    return -1;
  }

  // Declare the supported options.
//  po::options_description desc("lbl_dal");
//  desc.add_options()
//      ("help", "produce help message")
//      ("-i", po::value<std::string>(), "input")
//      ;
//
//  po::variables_map vm;
//  po::store(po::parse_command_line(argc, (const char *const *) argv, desc), vm);
//  po::notify(vm);
//
//  bool is_world_rank0 = p_ops->world_proc_rank == 0;
//  if (vm.count("help")) {
//    std::cout << desc << "\n";
//    return 1;
//  }
//
//  if (vm.count("-i")){
//    input_file = vm["-i"].as<std::string>();
//  }else {
//    if (is_world_rank0)
//      std::cout<<"ERROR: Input file not specified"<<std::endl;
//    return -1;
//  }

  return 0;
}
