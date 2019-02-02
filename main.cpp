#include <iostream>
#include <cxxopts.hpp>
#include <cmath>
#include "Constants.hpp"
#include "ParallelOps.hpp"
#include "ParallelFastaReader.hpp"
#include "Alphabet.hpp"
#include "Kmer.hpp"
#include "Utils.hpp"
#include "CombBLAS/CombBLAS.h"

/*! Namespace declarations */
using namespace combblas;

/*! Type definitions */
typedef KmerIntersect<int64_t, CommonKmers> KmerIntersectSR_t;

/*! Function signatures */
int parse_args(int argc, char **argv);

/*! Global variables */
std::shared_ptr<ParallelOps> p_ops;
std::string input_file;
int input_overlap;
int input_seq_count;
int klength;
int kstride;


int main(int argc, char **argv) {
  p_ops = ParallelOps::init(&argc, &argv);
  int ret = parse_args(argc, argv);
  if (ret < 0){
    p_ops->teardown_parallelism();
    return ret;
  }

  std::shared_ptr<FastaData> fd;
  ParallelFastaReader pfr;
  pfr.readFasta(input_file.c_str(), input_seq_count, input_overlap,
                p_ops->world_proc_rank, p_ops->world_procs_count, fd);

  Alphabet alph(Alphabet::PROTEIN);

  /* Find k-mers */
  std::shared_ptr<char> seq;
  int len, start_offset, end_offset_inclusive;
  //Note: Saliya - this can be parallelized using OpenMP
  std::vector<int64_t> lrow_ids, lcol_ids, lvals;
  for (int lseq_idx = 0; lseq_idx < fd->count(); ++lseq_idx){
    seq = fd->get_sequence(lseq_idx, false, len, start_offset,
        end_offset_inclusive);
//    std::printf("rank: %d len %d\n", p_ops->world_proc_rank, len);
    auto num_kmers = add_kmers(seq.get(), len, start_offset,
                               end_offset_inclusive, klength, kstride, alph,
                               lcol_ids, lvals);
    lrow_ids.insert(lrow_ids.end(), num_kmers, lseq_idx);
  }

  /*if (p_ops->world_proc_rank == 0) {
    for (int64_t &row_id : lrow_ids) {
      std::cout << row_id << " ";
    }
    std::cout << std::endl;
  }*/

//  assert(lrow_ids.size() == lcol_ids.size() == lvals.size());
  std::printf("Rank: %d lrow_ids %ld lcol_ids %ld lvals %ld\n",
              p_ops->world_proc_rank, lrow_ids.size(), lcol_ids.size(), lvals.size());


  auto grid_size = static_cast<int>(sqrt(p_ops->world_procs_count));
  std::shared_ptr<CommGrid> grid = std::make_shared<CommGrid>(
      MPI_COMM_WORLD, grid_size, grid_size);
  FullyDistVec<int64_t, int64_t> drows(lrow_ids, grid);
  FullyDistVec<int64_t, int64_t> dcols(lcol_ids, grid);
  FullyDistVec<int64_t, int64_t> dvals(lvals, grid);

  int64_t n_rows = input_seq_count;
  /*! Columns of the matrix are direct maps to kmers identified
   * by their |alphabet| base number. E.g. for proteins this is
   * base 20, so the direct map has to be 20^k in size.
   */
  auto n_cols = static_cast<int64_t>(pow(alph.size, klength));
  PSpMat<int64_t >::MPI_DCCols A(n_rows, n_cols, drows, dcols, dvals, false);

//  A.PrintInfo();

  auto At = A;
  At.Transpose();

  PSpMat<CommonKmers>::MPI_DCCols C =
      Mult_AnXBn_Synch<KmerIntersectSR_t,
          CommonKmers, PSpMat<CommonKmers>::DCCols>(A, At);



//  int nnz = C.getnnz();
//  if (p_ops->world_proc_rank == 0){
//    std::cout<<nnz<<std::endl;
//  }
//  C.PrintInfo();

  /*! Test multiplication */
  /* rows and cols in the result */

  int nrows = input_seq_count;
  int ncols = input_seq_count;
  int pr = static_cast<int>(sqrt(p_ops->world_procs_count));
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
  if (p_ops->world_proc_rank == 0) {
    for (auto colit = spSeq->begcol();
         colit != spSeq->endcol(); ++colit) // iterate over columns
    {
      std::cout<<"came here";
      int64_t lj = colit.colid(); // local numbering
      int64_t j = lj + localColStart;

      for (auto nzit = spSeq->begnz(colit);
           nzit < spSeq->endnz(colit); ++nzit) {

        int64_t li = nzit.rowid();
        int64_t i = li + localRowStart;
        std::cout<<"r:"<<li<<" c:"<<lj<<" v:"<<nzit.value()<<std::endl;
      }
    }
  }


  /* Test get fasta */
  /*int len, start_offset, end_offset_inclusive;
  std::shared_ptr<char> seq = fd->get_sequence(
      1, false, len, start_offset, end_offset_inclusive);

  int flag;
  if (p_ops->world_proc_rank > 0)
    MPI_Recv(&flag, 1, MPI_INT, p_ops->world_proc_rank-1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::cout<<"\nRank: "<<p_ops->world_proc_rank<<"\n-------------------------\n";
  std::cout.write(seq.get()+start_offset, (end_offset_inclusive - start_offset)+1);
  std::cout<< std::endl << start_offset << " " << end_offset_inclusive<<std::endl;
  if (p_ops->world_proc_rank < p_ops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, p_ops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
  }*/



  /* Test print */
  /*flag;
  if (p_ops->world_proc_rank > 0)
    MPI_Recv(&flag, 1, MPI_INT, p_ops->world_proc_rank-1, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  std::cout<<"\nRank: "<<p_ops->world_proc_rank<<"\n-------------------------\n";
  fd->print();
  if (p_ops->world_proc_rank < p_ops->world_procs_count - 1) {
    MPI_Send(&flag, 1, MPI_INT, p_ops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
  }*/

  p_ops->teardown_parallelism();
  return 0;
}

int parse_args(int argc, char **argv) {
  cxxopts::Options options("Distributed Aligner", "A distributed protein aligner");

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
       cxxopts::value<int>())
      ;

  auto result = options.parse(argc, argv);

  bool is_world_rank0 = p_ops->world_proc_rank == 0;
  if (result.count(CMD_OPTION_INPUT)){
    input_file = result[CMD_OPTION_INPUT].as<std::string>();
  }else {
    if (is_world_rank0) {
      std::cout << "ERROR: Input file not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_INPUT_SEQ_COUNT)){
    input_seq_count = result[CMD_OPTION_INPUT_SEQ_COUNT].as<int>();
  }else {
    if (is_world_rank0) {
      std::cout << "ERROR: Input sequence count not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_INPUT_OVERLAP)){
    input_overlap = result[CMD_OPTION_INPUT_OVERLAP].as<int>();
  }else {
    input_overlap = 10000;
  }

  if (result.count(CMD_OPTION_KMER_LENGTH)){
    klength = result[CMD_OPTION_KMER_LENGTH].as<int>();
  }else {
    if (is_world_rank0) {
      std::cout << "ERROR: Kmer length not specified" << std::endl;
    }
    return -1;
  }

  if (result.count(CMD_OPTION_KMER_STRIDE)){
    kstride = result[CMD_OPTION_KMER_STRIDE].as<int>();
  }else {
    if (is_world_rank0){
     kstride = 1;
    }
  }

  return 0;
}
