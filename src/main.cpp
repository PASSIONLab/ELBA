#include <iostream>
#include <cmath>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include "../include/Constants.hpp"
#include "../include/ParallelOps.hpp"
#include "../include/ParallelFastaReader.hpp"
#include "../include/Alphabet.hpp"
#include "../include/Kmer.hpp"
#include "../include/Utils.hpp"
#include "../include/DistributedPairwiseRunner.hpp"
#include "CombBLAS/CombBLAS.h"
#include "cxxopts.hpp"
#include "../include/pw/SeedExtendXdrop.hpp"
#include "seqan/score/score_matrix_data.h"
#include "../include/pw/OverlapFinder.hpp"
#include "../include/ScoreMat.hpp"
#include "../include/pw/FullAligner.hpp"
#include "../include/pw/BandedAligner.hpp"
#include <map>

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
uint64_t input_overlap;
uint64_t seq_count;
int xdrop;
int gap_open;
int gap_ext;
ushort klength;
ushort kstride;

/*! Parameters related to outputting k-mer overlaps */
bool write_overlaps = false;
std::string overlap_file;
bool add_substitue_kmers = false;
float subtitute_kmer_percentage = 100;

/*! Parameters related to outputting alignment info */
std::string align_file;

/*! Don't perform alignments if this flag is set */
bool no_align = false;

/*! Perform full alignment */
bool full_align = false;

/*! Perform xdrop alignment */
bool xdrop_align = false;

/*! Perform banded alignment */
bool banded_align = false;
int banded_half_width = 5;

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
std::string job_name = "";
std::string proc_log_file;
std::ofstream proc_log_stream;
int log_freq;

int main(int argc, char **argv) {
  parops = ParallelOps::init(&argc, &argv);
  int ret = parse_args(argc, argv);
  if (ret < 0) {
    parops->teardown_parallelism();
    return ret;
  }

  /*! bcast job id */
  // Sender
  int job_name_length = job_name.length();

  if(parops->world_proc_rank == 0) {
    MPI_Bcast(&job_name[0], job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
  } else {
    char* buf = new char[job_name_length];
    MPI_Bcast(buf, job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
    std::string s(buf);
    job_name = s;
    delete [] buf;
  }


  proc_log_file = job_name + "_rank_" + std::to_string(parops->world_proc_rank) + "_log.txt";
  proc_log_stream.open(proc_log_file);

  is_print_rank = (parops->world_proc_rank == 0);
  std::shared_ptr<TimePod> tp = std::make_shared<TimePod>();
  TraceUtils tu(is_print_rank);

  /*! Print start time information */
  tp->times["start_main"] = std::chrono::system_clock::now();
  std::time_t start_prog_time = std::chrono::system_clock::to_time_t(
    tp->times["start_main"]);
  print_str = "\nINFO: Program started on ";
  print_str.append(std::ctime(&start_prog_time));
  print_str.append("\nINFO: Job ID ").append(job_name).append("\n");
  pretty_print_config(print_str);
  tu.print_str(print_str);

  /*! Read and distribute fasta data */
  tp->times["start_main:newDFD()"] = std::chrono::system_clock::now();
  std::shared_ptr<DistributedFastaData> dfd
    = std::make_shared<DistributedFastaData>(
      input_file.c_str(), idx_map_file.c_str(), input_overlap,
      klength, parops, tp, tu);
  tp->times["end_main:newDFD()"] = std::chrono::system_clock::now();

#ifndef NDEBUG
  //  TraceUtils::print_fasta_data(fd, parops);
#endif

  if (dfd->global_count() != seq_count) {
    uint64_t final_seq_count = dfd->global_count();
    print_str = "\nINFO: Modfied sequence count\n";
    print_str.append("  Final sequence count: ")
      .append(std::to_string(final_seq_count))
      .append(" (").append(
        std::to_string((((seq_count - final_seq_count) * 100.0) / seq_count)))
      .append("% removed)");

    seq_count = dfd->global_count();
    print_str += "\n";
    tu.print_str(print_str);
  }

  /*! Find K-mers */
  Alphabet alph(alph_t);

  /*! Find k-mers */
  char *buff;
  ushort len;
  uint64_t start_offset, end_offset_inclusive;
  //TODO: Saliya - this can be parallelized using OpenMP
  uvec_64 lrow_ids, lcol_ids;
  uvec_16 lvals;
  uint64_t offset = dfd->global_start_idx();
  FastaData *lfd = dfd->lfd();
  pisa::Blosum62 bsm62;
  std::map<uint64_t, std::vector<uint64_t>*> kmer_to_subs_kmers;
  tp->times["start_main:loop_add_kmers()"] = std::chrono::system_clock::now();
  for (uint64_t lseq_idx = 0; lseq_idx < lfd->local_count(); ++lseq_idx) {
    buff = lfd->get_sequence(lseq_idx, len, start_offset, end_offset_inclusive);
    auto num_kmers = add_kmers(buff, len, start_offset, end_offset_inclusive,
                               klength, kstride, add_substitue_kmers,
                               subtitute_kmer_percentage,
                               kmer_to_subs_kmers, bsm62,
                               alph, lcol_ids, lvals, parops);
    lrow_ids.insert(lrow_ids.end(), num_kmers, lseq_idx + offset);
  }
  tp->times["end_main:loop_add_kmers()"] = std::chrono::system_clock::now();

  /*! Free kmer_to_subs_kmers vectors */
  for (std::pair<uint64_t, std::vector<uint64_t>*> p : kmer_to_subs_kmers){
    free(p.second);
  }


#ifndef NDEBUG
  {
    std::string title = "Local matrix info:";
    std::string msg = "lrow_ids row_size: " + std::to_string(lrow_ids.size())
                      + " lcol_ids row_size: " + std::to_string(lcol_ids.size())
                      + " lvals row_size: " + std::to_string(lvals.size());
    TraceUtils::print_msg(title, msg, parops);
  }
#endif

  assert(lrow_ids.size() == lcol_ids.size() && lcol_ids.size() == lvals.size());

  /*! Create distributed sparse matrix of sequence x kmers */
  FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, parops->grid);
  FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, parops->grid);
  /*! TODO - apparently there's a bug when setting a different element type,
   * so let's use uint64_t as element type for now */
  FullyDistVec<uint64_t, uint64_t> dvals(lvals, parops->grid);

  uint64_t n_rows = seq_count;
  /*! Columns of the matrix are direct maps to kmers identified
   * by their |alphabet| base number. E.g. for proteins this is
   * base 20, so the direct map has to be 20^k in size. */

  auto n_cols = static_cast<uint64_t>(pow(alph.size, klength));
  tp->times["start_main:spMatA()"] = std::chrono::system_clock::now();
  PSpMat<ushort>::MPI_DCCols A(n_rows, n_cols, drows, dcols, dvals, false);
  tp->times["end_main:spMatA()"] = std::chrono::system_clock::now();

//#ifndef NDEBUG
  A.PrintInfo();
//#endif

  auto At = A;
  tp->times["start_main:At()"] = tp->times["end_main:spMatA()"];
  At.Transpose();
  tp->times["end_main:At()"] = std::chrono::system_clock::now();

  tp->times["start_main:AxAt()"] = tp->times["end_main:At()"];
  PSpMat<CommonKmers>::MPI_DCCols C =
    Mult_AnXBn_Synch<KmerIntersectSR_t,
      CommonKmers, PSpMat<CommonKmers>::DCCols>(A, At);
  tu.print_str(
    "Overlaps after k-mer finding (nnz(C) - diagonal): "
    + std::to_string(C.getnnz() - seq_count)
    + "\nLoad imbalance: " + std::to_string(C.LoadImbalance()) + "\n");
  tp->times["end_main:AxAt()"] = std::chrono::system_clock::now();


#ifndef NDEBUG
  /*! Test multiplication */

  C.PrintInfo();

  // rows and cols in the result
  n_cols = n_rows = seq_count;
  int pr = parops->grid->GetGridRows();
  int pc = parops->grid->GetGridCols();

  int row_rank = parops->grid->GetRankInProcRow();
  int col_rank = parops->grid->GetRankInProcCol();
  uint64_t m_perproc = n_rows / pr;
  uint64_t n_perproc = n_cols / pc;
  PSpMat<CommonKmers>::DCCols *spSeq = C.seqptr(); // local submatrix
  uint64_t l_row_start = col_rank * m_perproc; // first row in this process
  uint64_t l_col_start = row_rank * n_perproc; // first col in this process

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


  std::cout << "\nRank: " << parops->world_proc_rank << "\n writing data\n";

  for (auto colit = spSeq->begcol();
       colit != spSeq->endcol(); ++colit) // iterate over columns
  {
    int64_t lj = colit.colid(); // local numbering
    int64_t j = lj + l_col_start;

    for (auto nzit = spSeq->begnz(colit);
         nzit < spSeq->endnz(colit); ++nzit) {

      int64_t li = nzit.rowid();
      int64_t i = li + l_row_start;
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
    MPI_Send(&flag, 1, MPI_INT, parops->world_proc_rank + 1, 99,
             MPI_COMM_WORLD);
  }

  MPI_Barrier(MPI_COMM_WORLD);
#endif

  /*! Wait until data distribution is complete */
  tp->times["start_main:dfd->wait()"] = std::chrono::system_clock::now();
  if (!dfd->is_ready()) {
    dfd->wait();
  }
  tp->times["end_main:dfd->wait()"] = std::chrono::system_clock::now();


  DistributedPairwiseRunner dpr(dfd, C, parops);
  if (!no_align) {
    tp->times["start_main:dpr->align()"] = std::chrono::system_clock::now();
    seqan::Blosum62 blosum62(gap_ext, gap_open);

    // TODO: SeqAn can't work with affine gaps for seed extension
    seqan::Blosum62 blosum62_simple(gap_open, gap_open);
    PairwiseFunction* pf = nullptr;
    uint64_t local_alignments = 1;
    if (xdrop_align) {
      pf = new SeedExtendXdrop (blosum62, blosum62_simple, klength, xdrop, seed_count);
      dpr.run(pf, align_file.c_str(), proc_log_stream, log_freq);
      local_alignments = static_cast<SeedExtendXdrop*>(pf)->alignments.size();
    } else if (full_align) {
      pf = new FullAligner(blosum62, blosum62_simple);
      dpr.run(pf, align_file.c_str(), proc_log_stream, log_freq);
      local_alignments = static_cast<FullAligner*>(pf)->alignments.size();
    } else if(banded_align){
      pf = new BandedAligner (blosum62, banded_half_width);
      dpr.run(pf, align_file.c_str(), proc_log_stream, log_freq);
      local_alignments = static_cast<BandedAligner*>(pf)->alignments.size();
    }

    tp->times["end_main:dpr->align()"] = std::chrono::system_clock::now();
    delete pf;

    uint64_t total_alignments = 0;
    MPI_Reduce(&local_alignments, &total_alignments, 1, MPI_UINT64_T, MPI_SUM, 0,
               MPI_COMM_WORLD);

    if (is_print_rank) {
      std::cout << "Final alignment (L+U-D) count: " << 2 * total_alignments
                << std::endl;
    }
  }

  tp->times["start_main:dpr->write_overlaps()"] = std::chrono::system_clock::now();
  if (write_overlaps){
    dpr.write_overlaps(overlap_file.c_str());
  }
  tp->times["end_main:dpr->write_overlaps()"] = std::chrono::system_clock::now();

  tp->times["end_main"] = std::chrono::system_clock::now();

  std::time_t end_prog_time = std::chrono::system_clock::to_time_t(
    tp->times["end_main"]);
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
  cxxopts::Options options("Distributed DistributedPairwiseRunner",
                           "A distributed protein aligner");

  options.add_options()
    (CMD_OPTION_INPUT, CMD_OPTION_DESCRIPTION_INPUT,
     cxxopts::value<std::string>())
    (CMD_OPTION_INPUT_SEQ_COUNT, CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT,
     cxxopts::value<int>())
    (CMD_OPTION_INPUT_OVERLAP, CMD_OPTION_DESCRIPTION_INPUT_OVERLAP,
     cxxopts::value<uint64_t>())
    (CMD_OPTION_SEED_COUNT, CMD_OPTION_DESCRIPTION_SEED_COUNT,
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
     cxxopts::value<float>())
    (CMD_OPTION_LOG_FREQ, CMD_OPTION_DESCRIPTION_LOG_FREQ,
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

  if (result.count(CMD_OPTION_GAP_OPEN)) {
    gap_open = result[CMD_OPTION_GAP_OPEN].as<int>();
  } else {
    gap_open = -11;
  }

  if (result.count(CMD_OPTION_GAP_EXT)) {
    gap_ext = result[CMD_OPTION_GAP_EXT].as<int>();
  } else {
    gap_ext = -2;
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
    if (tmp == "protein"){
      alph_t = Alphabet::PROTEIN;
    } else {
      if (is_world_rank0) {
        std::cout << "ERROR: Unsupported alphabet type " << tmp << std::endl;
      }
    }
  } else {
    alph_t = Alphabet::PROTEIN;
  }

  if (result.count(CMD_OPTION_SUBS)) {
    add_substitue_kmers = true;
    subtitute_kmer_percentage = result[CMD_OPTION_SUBS].as<float>();
  }

  boost::uuids::random_generator gen;
  boost::uuids::uuid id = gen();
  if (result.count(CMD_OPTION_JOB_NAME_PREFIX)) {
    std::string tmp = result[CMD_OPTION_JOB_NAME_PREFIX].as<std::string>();
    job_name = tmp + "_" +  boost::uuids::to_string(id);
  } else {
    job_name = boost::uuids::to_string(id);
  }

  if (result.count(CMD_OPTION_LOG_FREQ)) {
    log_freq = result[CMD_OPTION_LOG_FREQ].as<int>();
  }

  return 0;
}

auto bool_to_str = [](bool b){
  return b ? "True" : "False";
};

void pretty_print_config(std::string &append_to) {
  std::vector<std::string> params = {
    "Input file (-i)",
    "Original sequence count (-c)",
    "Kmer length (k)",
    "Kmer stride (s)",
    "Overlap in bytes (-O)",
    "Max seed count (--sc)",
    "Gap open penalty (-g)",
    "Gap extension penalty (-e)",
    "Overlap file (--of)",
    "Alignment file (--af)",
    "No align (--na)",
    "Full align (--fa)",
    "Xdrop align (--xa)",
    "Banded align (--ba)",
    "Index map (--idxmap)",
    "Alphabet (--alph)",
    "Use substitute kmers (--subs)"
  };



  std::vector<std::string> vals = {
    input_file,
    std::to_string(seq_count),
    std::to_string(klength),
    std::to_string(kstride),
    std::to_string(input_overlap),
    std::to_string(seed_count),
    std::to_string(gap_open),
    std::to_string(gap_ext),
    !overlap_file.empty() ? overlap_file : "None",
    !align_file.empty() ? align_file : "None",
    bool_to_str(no_align),
    bool_to_str(full_align),
    bool_to_str(xdrop_align) + (xdrop_align ? "| xdrop: " + std::to_string(xdrop) : ""),
    bool_to_str(banded_align) + (banded_align ? " | half band: " + std::to_string(banded_half_width) : ""),
    !idx_map_file.empty() ? idx_map_file : "None",
    std::to_string(alph_t),
    bool_to_str(add_substitue_kmers) + (add_substitue_kmers ? " | sub kmers: " + std::to_string(subtitute_kmer_percentage) : ""),
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
//        static_cast<ushort>(max_length + 1 - param.row_size()), prefix))
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
