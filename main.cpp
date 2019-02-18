#include <iostream>
#include <cmath>
#include "Constants.hpp"
#include "ParallelOps.hpp"
#include "ParallelFastaReader.hpp"
#include "Alphabet.hpp"
#include "Kmer.hpp"
#include "Utils.hpp"
//#include "DebugUtils.hpp"
#include "CombBLAS/CombBLAS.h"
#include "cxxopts.hpp"

/*! Namespace declarations */
using namespace combblas;

/*! Type definitions */
typedef KmerIntersect<ushort, CommonKmers> KmerIntersectSR_t;

/*! Function signatures */
int parse_args(int argc, char **argv);

void pretty_print_config(std::string &append_to);

std::string get_padding(ushort count, std::string prefix);

/*! Global variables */
std::shared_ptr<ParallelOps> parops;
std::string input_file;
ushort input_overlap;
uint64_t seq_count;
ushort klength;
ushort kstride;

bool is_print_rank = false;
std::string print_str;


int main(int argc, char **argv) {
  parops = ParallelOps::init(&argc, &argv);
  int ret = parse_args(argc, argv);
  if (ret < 0) {
    parops->teardown_parallelism();
    return ret;
  }

  is_print_rank = (parops->world_proc_rank == 0);

  /*! Print start time information */
  ticks_t start_prog = std::chrono::system_clock::now();
  std::time_t start_prog_time = std::chrono::system_clock::to_time_t(
    start_prog);
  print_str = "\nINFO: Program started on ";
  print_str.append(std::ctime(&start_prog_time));
  pretty_print_config(print_str);
  if (is_print_rank) {
    std::cout << print_str;
  }

  /*! Read and distribute fasta data */
  std::unique_ptr<DistributedFastaData> dfd
  = std::make_unique<DistributedFastaData>(
    input_file.c_str(), input_overlap, klength, parops);

#ifndef NDEBUG
//  DebugUtils::print_fasta_data(fd, parops);
#endif

  if (dfd->global_count() != seq_count) {
    uint64_t final_seq_count = dfd->global_count();
    if (is_print_rank) {
      print_str = "\nINFO: Modfied sequence count\n";
      print_str.append("  Final sequence count: ")
        .append(std::to_string(final_seq_count))
        .append(" (").append(
          std::to_string(((seq_count - final_seq_count) * 100 / seq_count)))
        .append("% removed)");

      seq_count = dfd->global_count();
      std::cout << print_str << std::endl;
    }
  }


  /*
  Alphabet alph(Alphabet::PROTEIN);

  *//* Find k-mers *//*
  std::shared_ptr<char> seq;
  ushort len;
  uint64_t start_offset, end_offset_inclusive;
  //TODO: Saliya - this can be parallelized using OpenMP
  uvec_64 lrow_ids, lcol_ids;
  uvec_16 lvals;
  uint64_t offset = fd->offset();
  for (uint64_t lseq_idx = 0; lseq_idx < fd->local_count(); ++lseq_idx){
    seq = fd->get_sequence(lseq_idx, false, len, start_offset,
        end_offset_inclusive);
//    std::printf("rank: %d len %d\n", parops->world_proc_rank, len);
    auto num_kmers = add_kmers(seq.get(), len, start_offset,
                               end_offset_inclusive, klength, kstride, alph,
                               lcol_ids, lvals);
    lrow_ids.insert(lrow_ids.end(), num_kmers, lseq_idx + offset);
  }


  *//*if (parops->world_proc_rank == 0) {
    for (int64_t &row_id : lrow_ids) {
      std::cout << row_id << " ";
    }
    std::cout << std::endl;
  }*//*

//  assert(lrow_ids.size() == lcol_ids.size() == lvals.size());
  std::printf("Rank: %d lrow_ids %ld lcol_ids %ld lvals %ld\n",
              parops->world_proc_rank, lrow_ids.size(), lcol_ids.size(), lvals.size());

  *//*! Write values to file to see why it segfaults *//*
*//*  std::ofstream f;
  std::string fname = "mat." + std::to_string(parops->world_proc_rank) + ".txt";
  f.open(fname);
  for (int i = 0; i < lrow_ids.size(); ++i){
    f << lrow_ids[i] << "," << lcol_ids[i] << "," << lvals[i] << std::endl;
  }
  f.close();*//*



  FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, grid);
  FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, grid);
  *//*! TODO - apparently there's a bug when setting a different element type,
   * so let's use uint64_t as element type for now
   *//*
  FullyDistVec<uint64_t, uint64_t> dvals(lvals, grid);

  uint64_t n_rows = seq_count;
  *//*! Columns of the matrix are direct maps to kmers identified
   * by their |alphabet| base number. E.g. for proteins this is
   * base 20, so the direct map has to be 20^k in size.
   *//*
  auto n_cols = static_cast<uint64_t>(pow(alph.size, klength));
  PSpMat<ushort>::MPI_DCCols A(n_rows, n_cols, drows, dcols, dvals, false);
  A.PrintInfo();

  auto At = A;
  At.Transpose();

  PSpMat<CommonKmers>::MPI_DCCols C =
      Mult_AnXBn_Synch<KmerIntersectSR_t,
          CommonKmers, PSpMat<CommonKmers>::DCCols>(A, At);

//  int nnz = C.getnnz();
//  if (parops->world_proc_rank == 0){
//    std::cout<<nnz<<std::endl;
//  }
  C.PrintInfo();

  *//*! Test multiplication *//*
  *//* rows and cols in the result *//*



  int64_t nrows = seq_count;
  int64_t ncols = seq_count;
  int pr = static_cast<int>(sqrt(parops->world_procs_count));
  int pc = pr;

  //Information about the matrix distribution
  //Assume that A is an nrow x ncol matrix
  //The local submatrix is an lnrow x lncol matrix
  int rowrank = grid->GetRankInProcRow();
  int colrank = grid->GetRankInProcCol();
  int64_t m_perproc = nrows / pr;
  int64_t n_perproc = ncols / pc;
  PSpMat<CommonKmers>::DCCols *spSeq = C.seqptr(); // local submatrix
  int64_t localRowStart = colrank * m_perproc; // first row in this process
  int64_t localColStart = rowrank * n_perproc; // first col in this process

//  std::cout<<std::endl<<spSeq->begcol().colid()<<" " <<spSeq->endcol().colid()<<std::endl;

  std::ofstream rf, cf, vf;
  int flag;
  if (parops->world_proc_rank > 0) {
    MPI_Recv(&flag, 1, MPI_INT, parops->world_proc_rank - 1,
             99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  rf.open("row_ids.txt", std::ios::app);
  cf.open("col_ids.txt", std::ios::app);
  vf.open("values.txt", std::ios::app);


  std::cout<<"\nRank: "<<parops->world_proc_rank<<"\n writing data\n";

  for (auto colit = spSeq->begcol();
       colit != spSeq->endcol(); ++colit) // iterate over columns
  {
    int64_t lj = colit.colid(); // local numbering
    int64_t j = lj + localColStart;

    for (auto nzit = spSeq->begnz(colit);
         nzit < spSeq->endnz(colit); ++nzit) {

      int64_t li = nzit.rowid();
      int64_t i = li + localRowStart;
      rf << i << ",";
      cf << j << ",";
      vf << nzit.value().count << ",";
//      std::cout<<"r:"<<li<<" c:"<<lj<<" v:"<<nzit.value()<<std::endl;
    }
  }

  rf.close();
  cf.close();
  vf.close();

  if (parops->world_proc_rank < parops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, parops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);*/

  /* Test get fasta */
  /*int len, start_offset, end_offset_inclusive;
  std::shared_ptr<char> seq = fd->get_sequence(
      1, false, len, start_offset, end_offset_inclusive);

  int flag;
  if (parops->world_proc_rank > 0)
    MPI_Recv(&flag, 1, MPI_INT, parops->world_proc_rank-1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::cout<<"\nRank: "<<parops->world_proc_rank<<"\n-------------------------\n";
  std::cout.write(seq.get()+start_offset, (end_offset_inclusive - start_offset)+1);
  std::cout<< std::endl << start_offset << " " << end_offset_inclusive<<std::endl;
  if (parops->world_proc_rank < parops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, parops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
  }



  /* Test print */
  /*flag;
  if (parops->world_proc_rank > 0)
    MPI_Recv(&flag, 1, MPI_INT, parops->world_proc_rank-1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::cout<<"\nRank: "<<parops->world_proc_rank<<"\n-------------------------\n";
  fd->print();
  if (parops->world_proc_rank < parops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, parops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
  }*/

  parops->teardown_parallelism();
  return 0;
}

int parse_args(int argc, char **argv) {
  cxxopts::Options options("Distributed Aligner",
                           "A distributed protein aligner");

  options.add_options()
    (CMD_OPTION_INPUT, CMD_OPTION_DESCRIPTION_INPUT,
     cxxopts::value<std::string>())
    (CMD_OPTION_INPUT_SEQ_COUNT, CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT,
     cxxopts::value<int>())
    (CMD_OPTION_INPUT_OVERLAP, CMD_OPTION_DESCRIPTION_INPUT_OVERLAP,
     cxxopts::value<int>())
    (CMD_OPTION_KMER_LENGTH, CMD_OPTION_DESCRIPTION_KMER_LENGTH,
     cxxopts::value<int>())
    (CMD_OPTION_KMER_STRIDE, CMD_OPTION_DESCRIPTION_KMER_STRID,
     cxxopts::value<int>());

  auto result = options.parse(argc, argv);

  bool is_world_rank0 = parops->world_proc_rank == 0;
  if (result.count(CMD_OPTION_INPUT)) {
    input_file = result[CMD_OPTION_INPUT].as<std::string>();
  } else {
    if (is_world_rank0) {
      std::cout << "ERROR: Input file not specified" << std::endl;
    }
    return -1;
  }

  // TODO - fix option parsing
  if (result.count(CMD_OPTION_INPUT_SEQ_COUNT)) {
    seq_count = result[CMD_OPTION_INPUT_SEQ_COUNT].as<int>();
  } else {
    if (is_world_rank0) {
      std::cout << "ERROR: Input sequence count not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_INPUT_OVERLAP)) {
    input_overlap = result[CMD_OPTION_INPUT_OVERLAP].as<int>();
  } else {
    input_overlap = 10000;
  }

  if (result.count(CMD_OPTION_KMER_LENGTH)) {
    klength = result[CMD_OPTION_KMER_LENGTH].as<int>();
  } else {
    if (is_world_rank0) {
      std::cout << "ERROR: Kmer length not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_KMER_STRIDE)) {
    kstride = result[CMD_OPTION_KMER_STRIDE].as<int>();
  } else {
    if (is_world_rank0) {
      kstride = 1;
    }
  }

  return 0;
}

void pretty_print_config(std::string &append_to) {
  std::vector<std::string> params = {
    "Input file (-i)",
    "Original sequence count (-c)",
    "Kmer length (k)",
    "Kmer stride (s)",
    "Overlap in bytes (-O)"
  };

  std::vector<std::string> vals = {
    input_file, std::to_string(seq_count),
    std::to_string(klength), std::to_string(kstride),
    std::to_string(input_overlap)
  };

  ushort max_length = 0;
  for (const auto &param : params) {
    if (param.size() > max_length) {
      max_length = static_cast<ushort>(param.size());
    }
  }

  std::string prefix = "  ";
  append_to.append("Parameters...\n");
  for (ushort i = 0; i < params.size(); ++i) {
    std::string param = params[i];
    append_to.append(prefix).append(param).append(": ")
      .append(get_padding(
        static_cast<ushort>(max_length - param.size()), ""))
      .append(vals[i]).append("\n");

//    append_to.append(get_padding(
//        static_cast<ushort>(max_length + 1 - param.size()), prefix))
//        .append(param).append(": ")
//        .append(vals[i]).append("\n");
  }
}

std::string get_padding(ushort count, std::string prefix) {
  std::string pad = std::move(prefix);
  for (ushort i = 0; i < count; ++i) {
    pad += " ";
  }
  return pad;
}
