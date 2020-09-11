#include <iostream>
#include <cmath>
#include "../include/Constants.hpp"
#include "../include/ParallelOps.hpp"
#include "../include/ParallelFastaReader.hpp"
#include "../include/Alphabet.hpp"
#include "../include/DistributedPairwiseRunner.hpp"
#include "../include/cxxopts.hpp"
#include "../include/pw/SeedExtendXdrop.hpp"
#include "../include/pw/OverlapFinder.hpp"
#include "../include/pw/FullAligner.hpp"
#include "../include/pw/BandedAligner.hpp"
#include "../include/kmer/KmerOps.hpp"
#include "../include/kmer/KmerIntersectSR.hpp"
#include "../include/TransitiveReductionSR.hpp"
#include "../include/Utils.hpp"

#include "seqan/score/score_matrix_data.h"

#include <map>
#include <fstream>

#define TWOSEED

/*! Namespace declarations */
using namespace combblas;

/*! Type definitions */
typedef dibella::KmerIntersect<PosInRead, dibella::CommonKmers> KmerIntersectSR_t;
typedef MinPlusBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> MinPlusSR_t;
typedef ReduceMBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> ReduceMSR_t;
typedef Bind2ndBiSRing <dibella::CommonKmers, dibella::CommonKmers, dibella::CommonKmers> Bind2ndSR_t;

/*! Function signatures */
int parse_args(int argc, char **argv);

void pretty_print_config(std::string &append_to);

std::string get_padding(ushort count, std::string prefix);

/*! Global variables */
std::shared_ptr<ParallelOps> parops;
std::string input_file;
uint64_t input_overlap;
uint64_t seq_count;
int xdrop;
int match;
int  mismatch_sc;
int gap_open;
int gap_ext;
ushort klength;
ushort kstride;

/*! Parameters related to outputting k-mer overlaps */
bool write_overlaps = false;
std::string overlap_file;
bool add_substitue_kmers = false;
int subk_count = 0;

/*! Parameters related to outputting alignment info */
std::string align_file;
int afreq;

/*! Don't perform alignments if this flag is set */
bool no_align = false;

/*! Perform full alignment */
bool full_align = false;

/*! Perform xdrop alignment */
bool xdrop_align = false;

/*! Perform banded alignment */
bool banded_align = false;
int  banded_half_width = 5;

/*! File path to output global sequence index to original global sequence
 * index mapping */
std::string idx_map_file;

/*! Alphabet to use. */
Alphabet::type alph_t;

bool is_print_rank = false;
std::string print_str;

/*! Maximum number of common k-mers to keep */
int seed_count = 2;

/*! Logging information */
std::string job_name = "dibella";
std::string proc_log_file;
std::ofstream proc_log_stream;
int log_freq;

/*! Common k-mer threshold */
int ckthr = 1;

/*! Score threshold */
bool aln_score_thr = false; // GGGG: Currently not used but I'm gonna leave it here for now in case it comes handy

int main(int argc, char **argv)
{
  parops = ParallelOps::init(&argc, &argv);
  int ret = parse_args(argc, argv);

  if (ret < 0)
  {
    parops->teardown_parallelism();
    return ret;
  }

  /*! bcast job id */
  // Sender
  int job_name_length = job_name.length();

  if(parops->world_proc_rank == 0)
  {
    MPI_Bcast(&job_name[0], job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
  }
  else
  {
    char* buf = new char[job_name_length];
    MPI_Bcast(buf, job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
    std::string s(buf);
    job_name = s;
    delete [] buf;
  }

    int nthreads = 1;
#ifdef THREADED
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
#endif
    
    int nprocs, myrank;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  
    if(myrank == 0)
    {
        std::cout << "Process Grid (p x p x t): " << sqrt(nprocs) << " x " << sqrt(nprocs) << " x " << nthreads << std::endl;
    }

  proc_log_file = job_name + "_rank_" + std::to_string(parops->world_proc_rank) + "_log.txt";
  proc_log_stream.open(proc_log_file);

  is_print_rank = (parops->world_proc_rank == 0);
  std::shared_ptr<TimePod> tp = std::make_shared<TimePod>();
  TraceUtils tu(is_print_rank);

  /*! Print start time information */
  tp->times["StartMain"] = std::chrono::system_clock::now();
  std::time_t start_prog_time = std::chrono::system_clock::to_time_t(
  tp->times["StartMain"]);
  
  print_str = "\nINFO: Program started on ";
  print_str.append(std::ctime(&start_prog_time));
  print_str.append("\nINFO: Job ID ").append(job_name).append("\n");
  pretty_print_config(print_str);
  tu.print_str(print_str);

  /*! Read and distribute fasta data */
  tp->times["StartMain:newDFD()"] = std::chrono::system_clock::now();

  std::shared_ptr<DistributedFastaData> dfd = std::make_shared<DistributedFastaData>(
      input_file.c_str(), idx_map_file.c_str(), input_overlap,
      klength, parops, tp, tu);
      
  tp->times["EndMain:newDFD()"] = std::chrono::system_clock::now();

#ifndef NDEBUG
  //  TraceUtils::print_fasta_data(fd, parops);
#endif

  if (dfd->global_count() != seq_count)
  {
    uint64_t final_seq_count = dfd->global_count();
    print_str = "\nINFO: Modified sequence count\n";
    print_str.append("  Final sequence count: ")
      .append(std::to_string(final_seq_count))
      .append(" (").append(
        std::to_string((((seq_count - final_seq_count) * 100.0) / seq_count)))
      .append("% removed)");

    seq_count = dfd->global_count();
    print_str += "\n";
    tu.print_str(print_str);
  }

  /*! Create alphabet */
  Alphabet alph(alph_t);

  tp->times["StartMain:GenerateA()"] = std::chrono::system_clock::now();
  PSpMat<PosInRead>::MPI_DCCols A =
      dibella::KmerOps::GenerateA(
          seq_count, dfd, klength, kstride,
          alph, parops, tp);

  tu.print_str("Matrix A: ");
  tu.print_str("\nLoad imbalance: " + std::to_string(A.LoadImbalance()) + "\n");

  tp->times["EndMain:GenerateA()"] = std::chrono::system_clock::now();

  A.PrintInfo();

  auto At = A;
  tp->times["StartMain:At()"] = tp->times["EndMain:GenerateA()"];
  At.Transpose();
  tu.print_str("Matrix At: ");
  At.PrintInfo();
  tp->times["EndMain:At()"] = std::chrono::system_clock::now();

  tp->times["StartMain:AAt()"] = std::chrono::system_clock::now();
  proc_log_stream << "INFO: Rank: " << parops->world_proc_rank << " starting AAt" << std::endl;

  // GGGG: there's some bug in the vector version (new one stack error)
  PSpMat<dibella::CommonKmers>::MPI_DCCols B = Mult_AnXBn_DoubleBuff<KmerIntersectSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(A, At);  

  proc_log_stream << "INFO: Rank: " << parops->world_proc_rank << " done AAt" << std::endl;
  tu.print_str(
      "Matrix AAt: Overlaps after k-mer finding (nnz(C) - diagonal): "
      + std::to_string(B.getnnz() - seq_count)
      + "\nLoad imbalance: " + std::to_string(B.LoadImbalance()) + "\n");
  tp->times["EndMain:AAt()"] = std::chrono::system_clock::now();

  tu.print_str("Matrix B, i.e AAt: ");
  B.PrintInfo();

  /*! Wait until data distribution is complete */
  tp->times["StartMain:DfdWait()"] = std::chrono::system_clock::now();
  if (!dfd->is_ready())
  {
    dfd->wait();
  }
  tp->times["EndMain:DfdWait()"] = std::chrono::system_clock::now();

  uint64_t n_rows, n_cols;
  n_rows = n_cols = dfd->global_count();
  int gr_rows = parops->grid->GetGridRows();
  int gr_cols = parops->grid->GetGridCols();

  int gr_col_idx = parops->grid->GetRankInProcRow();
  int gr_row_idx = parops->grid->GetRankInProcCol();

  uint64_t avg_rows_in_grid = n_rows / gr_rows;
  uint64_t avg_cols_in_grid = n_cols / gr_cols;
  uint64_t row_offset = gr_row_idx * avg_rows_in_grid;  // first row in this process
  uint64_t col_offset = gr_col_idx * avg_cols_in_grid;	// first col in this process

  DistributedPairwiseRunner dpr(dfd, B.seqptr(), &B, afreq, row_offset, col_offset, parops);

  if (!no_align)
  {

    tp->times["StartMain:DprAlign()"] = std::chrono::system_clock::now();
    ScoringScheme scoring_scheme(match, mismatch_sc, gap_ext);

    align_file += "_rank_" + std::to_string(parops->world_proc_rank) + ".txt";

    PairwiseFunction* pf = nullptr;
    uint64_t local_alignments = 1;

    if (xdrop_align)
    {
      pf = new SeedExtendXdrop (scoring_scheme, klength, xdrop, seed_count);	    
      dpr.run_batch(pf, align_file, proc_log_stream, log_freq, ckthr, aln_score_thr, tu, klength);
	    local_alignments = static_cast<SeedExtendXdrop*>(pf)->nalignments;
    }
    else if (full_align)
    {
      pf = new FullAligner(scoring_scheme);
      dpr.run_batch(pf, align_file, proc_log_stream, log_freq, ckthr, aln_score_thr, tu, klength);
	    local_alignments = static_cast<FullAligner*>(pf)->nalignments;
    }
    else if(banded_align)
    {
      pf = new BandedAligner (scoring_scheme, banded_half_width);
      dpr.run_batch(pf, align_file, proc_log_stream, log_freq, ckthr, aln_score_thr, tu, klength);
	    local_alignments = static_cast<BandedAligner*>(pf)->nalignments;
    }

    tp->times["EndMain:DprAlign()"] = std::chrono::system_clock::now();
    delete pf;
	
    uint64_t total_alignments = 0;
    MPI_Reduce(&local_alignments, &total_alignments, 1, MPI_UINT64_T, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (is_print_rank)
    {
      std::cout << "Final alignment (L+U-D) count: " << 2 * total_alignments
                << std::endl;
    }
  }

  tp->times["StartMain:DprWriteOverlaps()"] = std::chrono::system_clock::now();
  if (write_overlaps)
  {
    dpr.write_overlaps(overlap_file.c_str());
  }
  tp->times["EndMain:DprWriteOverlaps()"] = std::chrono::system_clock::now();

  tp->times["StartMain:TransitiveReduction()"] = std::chrono::system_clock::now();
  bool transitive_reduction = true;
  if (transitive_reduction)
  {
    /* Randomly permute for load balance */
    // int perm = false; // GGGG: TODO make this input parameter
    // if(perm == 1)
    // {
    //     SpParHelper::Print("Performing random permutation of matrix\n");

    //     FullyDistVec<int64_t, dibella::CommonKmers> prow(A.getcommgrid());
    //     FullyDistVec<int64_t, dibella::CommonKmers> pcol(A.getcommgrid());

    //     // GGGG: from line 931 of FullyDistVect.cpp:
    //     // ABAB: In its current form, unless LengthUntil returns NT
    //     // this won't work reliably for anything other than integers
    //     // GGGG: TODO look into this
    //
    //     prow.iota(A.getnrow(), id);
    //     pcol.iota(A.getncol(), id);

    //     prow.RandPerm();
    //     pcol.RandPerm();
    //     (A)(prow, pcol, true);

    //     SpParHelper::Print("Performed random permutation of matrix\n");
    // }

    uint nnz, prev;

    /* Gonna iterate on B until there are no more transitive edges to remove */
    do
    {
      prev = B.getnnz();

      /* Find two-hops neighbors
       * C = B^2
       */
      PSpMat<dibella::CommonKmers>::MPI_DCCols F = B;
      PSpMat<dibella::CommonKmers>::MPI_DCCols C = Mult_AnXBn_DoubleBuff<MinPlusSR_t, dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, F);
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix C = B^2: ");
      C.PrintInfo();
  #endif
    
      FullyDistVec<int64_t, dibella::CommonKmers> vA(B.getcommgrid());

      dibella::CommonKmers id; // Double check id() here
      vA = B.Reduce(Row, ReduceMSR_t(), id);
      vA.Apply(PlusFBiSRing<dibella::CommonKmers, dibella::CommonKmers>());

      F.DimApply(Row, vA, Bind2ndSR_t());
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix F = B + FUZZ: ");
      F.PrintInfo();
  #endif
      /* Find transitive edges that can be removed
       * I = F >= C 
       */
      bool isLogicalNot = false;
      PSpMat<bool>::MPI_DCCols I = EWiseApply<bool, PSpMat<bool>::DCCols>(F, C, GreaterBinaryOp<dibella::CommonKmers, dibella::CommonKmers>(), isLogicalNot, id);

      /* Prune potential zero-valued nonzeros */
      I.Prune(ZeroUnaryOp<bool>(), true);
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix I = F >= B: ");
      I.PrintInfo();
  #endif

      /* Remove transitive edges
       * B = B .* not(I)
       */ 
      isLogicalNot = true;
      B = EWiseApply<dibella::CommonKmers, PSpMat<dibella::CommonKmers>::DCCols>(B, I, EWiseMulOp<dibella::CommonKmers, bool>(), isLogicalNot, true);

      /* Prune zero-valued overhang */
      B.Prune(ZeroOverhangSR<dibella::CommonKmers>(), true);
  #ifdef DIBELLA_DEBUG
      tu.print_str("Matrix B = B .* not(I): ");
      B.PrintInfo();
  #endif
      nnz = B.getnnz();
       
    } while (nnz != prev);
  }
  tu.print_str("Matrix B, i.e AAt after transitive reduction: ");
  B.PrintInfo();
  tp->times["EndMain:TransitiveReduction()"] = std::chrono::system_clock::now();

  tp->times["EndMain"] = std::chrono::system_clock::now();

  std::time_t end_prog_time = std::chrono::system_clock::to_time_t(
    tp->times["EndMain"]);
  print_str = "INFO: Program ended on ";
  print_str.append(std::ctime(&end_prog_time));
  tu.print_str(print_str);
  tu.print_str(tp->to_string());

  proc_log_stream.close();

  MPI_Barrier(MPI_COMM_WORLD);
  parops->teardown_parallelism();
  return 0;
}

int parse_args(int argc, char **argv) {
  cxxopts::Options options("diBELLA",
                           "Distributed Long Read to Long Read Alignment");

  options.add_options()
    (CMD_OPTION_INPUT, CMD_OPTION_DESCRIPTION_INPUT,
     cxxopts::value<std::string>())
    (CMD_OPTION_INPUT_SEQ_COUNT, CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT,
     cxxopts::value<int>())
    (CMD_OPTION_INPUT_OVERLAP, CMD_OPTION_DESCRIPTION_INPUT_OVERLAP,
     cxxopts::value<uint64_t>())
    (CMD_OPTION_SEED_COUNT, CMD_OPTION_DESCRIPTION_SEED_COUNT,
     cxxopts::value<int>())
    (CMD_OPTION_MATCH, CMD_OPTION_DESCRIPTION_MATCH,
    cxxopts::value<int>())
    (CMD_OPTION_MISMATCH, CMD_OPTION_DESCRIPTION_MISMATCH,
    cxxopts::value<int>())
    (CMD_OPTION_GAP_OPEN, CMD_OPTION_DESCRIPTION_GAP_OPEN,
     cxxopts::value<int>())
    (CMD_OPTION_GAP_EXT, CMD_OPTION_DESCRIPTION_GAP_EXT,
     cxxopts::value<int>())
    (CMD_OPTION_KMER_LENGTH, CMD_OPTION_DESCRIPTION_KMER_LENGTH,
     cxxopts::value<int>())
    (CMD_OPTION_KMER_STRIDE, CMD_OPTION_DESCRIPTION_KMER_STRID,
     cxxopts::value<int>())
    (CMD_OPTION_OVERLAP_FILE, CMD_OPTION_DESCRIPTION_OVERLAP_FILE,
     cxxopts::value<std::string>())
    (CMD_OPTION_ALIGN_FILE, CMD_OPTION_DESCRIPTION_ALIGN_FILE,
     cxxopts::value<std::string>())
    (CMD_OPTION_NO_ALIGN, CMD_OPTION_DESCRIPTION_NO_ALIGN)
    (CMD_OPTION_FULL_ALIGN, CMD_OPTION_DESCRIPTION_FULL_ALIGN)
    (CMD_OPTION_XDROP_ALIGN, CMD_OPTION_DESCRIPTION_XDROP_ALIGN,
     cxxopts::value<int>())
    (CMD_OPTION_BANDED_ALIGN, CMD_OPTION_DESCRIPTION_BANDED_ALIGN,
     cxxopts::value<int>())
    (CMD_OPTION_IDX_MAP, CMD_OPTION_DESCRIPTION_IDX_MAP,
     cxxopts::value<std::string>())
    (CMD_OPTION_ALPH, CMD_OPTION_DESCRIPTION_ALPH,
     cxxopts::value<std::string>())
    (CMD_OPTION_JOB_NAME_PREFIX, CMD_OPTION_DESCRIPTION_JOB_NAME_PREFIX,
     cxxopts::value<std::string>())
    (CMD_OPTION_SUBS, CMD_OPTION_DESCRIPTION_SUBS,
     cxxopts::value<int>())
    (CMD_OPTION_LOG_FREQ, CMD_OPTION_DESCRIPTION_LOG_FREQ,
     cxxopts::value<int>())
    (CMD_OPTION_AF_FREQ, CMD_OPTION_DESCRIPTION_AF_FREQ,
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

    if (result.count(CMD_OPTION_IDX_MAP)) {
        idx_map_file = result[CMD_OPTION_IDX_MAP].as<std::string>();
    } else {
        if (is_world_rank0) {
            std::cout << "ERROR: Index map file not specified" << std::endl;
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
    input_overlap = result[CMD_OPTION_INPUT_OVERLAP].as<uint64_t>();
  } else {
    input_overlap = 10000;
  }

  if (result.count(CMD_OPTION_SEED_COUNT)) {
    seed_count = result[CMD_OPTION_SEED_COUNT].as<int>();
  } else {
    seed_count = 2;
  }

  if (result.count(CMD_OPTION_MATCH)) {
    match = result[CMD_OPTION_MATCH].as<int>();
  } else {
    match = 1;
  }

  if (result.count(CMD_OPTION_MISMATCH)) {
     mismatch_sc = result[CMD_OPTION_MISMATCH].as<int>();
  } else {
     mismatch_sc = -1;
  }

  if (result.count(CMD_OPTION_GAP_OPEN)) {
    gap_open = result[CMD_OPTION_GAP_OPEN].as<int>();
  } else {
    gap_open = 0;
  }

  if (result.count(CMD_OPTION_GAP_EXT)) {
    gap_ext = result[CMD_OPTION_GAP_EXT].as<int>();
  } else {
    gap_ext = -1;
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

  if (result.count(CMD_OPTION_OVERLAP_FILE)) {
    write_overlaps = true;
    overlap_file = result[CMD_OPTION_OVERLAP_FILE].as<std::string>();
  }

  if (result.count(CMD_OPTION_ALIGN_FILE)) {
    align_file = result[CMD_OPTION_ALIGN_FILE].as<std::string>();
  }

  if (result.count(CMD_OPTION_NO_ALIGN)) {
    no_align = true;
  }

  if (result.count(CMD_OPTION_FULL_ALIGN)) {
    full_align = true;
  }

  if (result.count(CMD_OPTION_BANDED_ALIGN)) {
    banded_align = true;
    banded_half_width = result[CMD_OPTION_BANDED_ALIGN].as<int>();
  }

  if (result.count(CMD_OPTION_XDROP_ALIGN)) {
    xdrop_align = true;
    xdrop = result[CMD_OPTION_XDROP_ALIGN].as<int>();
  }

  if (result.count(CMD_OPTION_ALPH)) {
    std::string tmp = result[CMD_OPTION_ALPH].as<std::string>();
    if (tmp == "dna"){
      alph_t = Alphabet::DNA;
    } else {
      if (is_world_rank0) {
        std::cout << "ERROR: Unsupported alphabet type " << tmp << std::endl;
      }
    }
  } else {
    alph_t = Alphabet::DNA;
  }

  if (result.count(CMD_OPTION_LOG_FREQ)) {
    log_freq = result[CMD_OPTION_LOG_FREQ].as<int>();
  }

  if (result.count(CMD_OPTION_AF_FREQ)) {
    afreq = result[CMD_OPTION_AF_FREQ].as<int>();
  }

  return 0;
}

auto bool_to_str = [](bool b){
  return b ? "True" : "False";
};

/*! GGGG: add input match/ mismatch_sc */
void pretty_print_config(std::string &append_to) {
  std::vector<std::string> params = {
    "Input file (-i)",
    "Original sequence count (-c)",
    "Kmer length (k)",
    "Kmer stride (s)",
    "Overlap in bytes (-O)",
    "Max seed count (--sc)",
    "Base match score (--ma)",
    "Base mismatch score (--mi)",
    "Gap open penalty (-g)",
    "Gap extension penalty (-e)",
    "Overlap file (--of)",
    "Alignment file (--af)",
    "Alignment write frequency (--afreq)",
    "No align (--na)",
    "Full align (--fa)",
    "Xdrop align (--xa)",
    "Banded align (--ba)",
    "Index map (--idxmap)",
    "Alphabet (--alph)"
  };

  std::vector<std::string> vals = {
    input_file,
    std::to_string(seq_count),
    std::to_string(klength),
    std::to_string(kstride),
    std::to_string(input_overlap),
    std::to_string(seed_count),
    std::to_string(match),
    std::to_string(mismatch_sc),
    std::to_string(gap_open),
    std::to_string(gap_ext),
    !overlap_file.empty() ? overlap_file  : "None",
    !align_file.empty()   ? align_file    : "None",
    !align_file.empty()   ? std::to_string(afreq) : "None",
    bool_to_str(no_align),
    bool_to_str(full_align),
    bool_to_str(xdrop_align)  + (xdrop_align  ? " | xdrop: " + std::to_string(xdrop) : ""),
    bool_to_str(banded_align) + (banded_align ? " | half band: " + std::to_string(banded_half_width) : ""),
    !idx_map_file.empty() ? idx_map_file : "None",
    std::to_string(alph_t)
  };

  ushort max_length = 0;
  for (const auto &param : params)
  {
    if (param.size() > max_length)
    {
      max_length = static_cast<ushort>(param.size());
    }
  }

  std::string prefix = "  ";
  append_to.append("Parameters:\n");

  for (ushort i = 0; i < params.size(); ++i)
  {
    std::string param = params[i];
    append_to.append(prefix).append(param).append(": ")
      .append(get_padding(
        static_cast<ushort>(max_length - param.size()), ""))
      .append(vals[i]).append("\n");
  }
}

std::string get_padding(ushort count, std::string prefix) {
  std::string pad = std::move(prefix);
  for (ushort i = 0; i < count; ++i) {
    pad += " ";
  }
  return pad;
}
