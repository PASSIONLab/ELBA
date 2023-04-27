/*
PASTIS Copyright (c) 2020, The Regents of the University of California, through Lawrence Berkeley National Laboratory
(subject to receipt of any required approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software, please contact Berkeley Lab's Intellectual
Property Office at IPO@lbl.gov.

NOTICE. This Software was developed under funding from the U.S. Department of Energy and the U.S. Government consequently
retains certain rights. As such, the U.S. Government has been granted for itself and others acting on its behalf a paid-up,
nonexclusive, irrevocable, worldwide license in the Software to reproduce, distribute copies to the public, prepare
derivative works, and perform publicly and display publicly, and to permit others to do so.
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdint>
#include <cassert>
#include <unistd.h>
#include <mpi.h>
#include "common.h"
#include "compiletime.h"
#include "Logger.hpp"
#include "FastaIndex.hpp"
#include "FastaData.hpp"
#include "Kmer/Kmer.hpp"

int returncode;
std::string fasta_fname;

/*
 * MPI communicator info.
 */
MPI_Comm comm;
int myrank;
int nprocs;

/*
 * X-Drop alignment parameters.
 */
int mat = 1;           /* match score */
int mis = -1;          /* mismatch score */
int gap = -1;          /* gap score */
int xdrop_cutoff = 15; /* x-drop cutoff score */

constexpr int root = 0; /* root process rank */

/*
 * Timer variable(s).
 */
double elapsed, mintime, maxtime, avgtime;

/*
 * Command line parser.
 */
int parse_cli(int argc, char *argv[]);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    {
        /*
         * CombBLAS 2D processor grid uses MPI_COMM_WORLD
         * as its world communicator (comm).
         */
        Grid commgrid;

        comm = MPI_COMM_WORLD;
        commgrid.reset(new CommGrid(comm, 0, 0));
        myrank = commgrid->GetRank();
        nprocs = commgrid->GetSize();

        /*
         * Parse command-line from root process and broadcast.
         */
        if (parse_cli(argc, argv) != 0)
            goto err;

        MPI_Barrier(comm);
        elapsed = -MPI_Wtime();

        /*
         * Start pipeline
         */

        FIndex index(new FastaIndex(fasta_fname, commgrid));
        FastaData lfd(index);
        lfd.log();

        /*
         * Finish pipeline
         */

        elapsed += MPI_Wtime();
        MPI_Reduce(&elapsed, &maxtime, 1, MPI_DOUBLE, MPI_MAX, root, comm);
        MPI_Reduce(&elapsed, &avgtime, 1, MPI_DOUBLE, MPI_SUM, root, comm);

        if (myrank == root)
        {
            std::cout << "Total time (user secs): " << std::fixed << std::setprecision(3) << elapsed << "\n"
                      << "Total work (proc secs): " << std::fixed << std::setprecision(3) << avgtime << "\n"
                      << "Total cost (proc secs): " << std::fixed << std::setprecision(3) << (elapsed*nprocs) << std::endl;
        }

        MPI_Barrier(comm);

        goto done;
    }

err: returncode = -1;
done: MPI_Finalize();
      return returncode;
}

void usage(char const *prg)
{
    std::cerr << "Usage: " << prg << " [options] <reads.fa>\n"
              << "Options: -x INT   x-drop alignment threshold [" <<  xdrop_cutoff << "]\n"
              << "         -A INT   matching score ["             <<  mat          << "]\n"
              << "         -B INT   mismatch penalty ["           << -mis          << "]\n"
              << "         -G INT   gap penalty ["                << -gap          << "]\n"
              << "         -h       help message"
              << std::endl;
}

int parse_cli(int argc, char *argv[])
{
    int params[4] = {mat, mis, gap, xdrop_cutoff};
    int show_help = 0, fasta_provided = 1;

    if (myrank == root)
    {
        int c;

        while ((c = getopt(argc, argv, "x:A:B:G:h")) >= 0)
        {
            if      (c == 'A') params[0] =  atoi(optarg);
            else if (c == 'B') params[1] = -atoi(optarg);
            else if (c == 'G') params[2] = -atoi(optarg);
            else if (c == 'x') params[3] =  atoi(optarg);
            else if (c == 'h') show_help = 1;
        }
    }

    MPI_BCAST(params, 4, MPI_INT, root, comm);

    mat          = params[0];
    mis          = params[1];
    gap          = params[2];
    xdrop_cutoff = params[3];

    if (myrank == root && show_help)
        usage(argv[0]);

    MPI_BCAST(&show_help, 1, MPI_INT, root, comm);
    if (show_help) return -1;

    if (myrank == root && optind >= argc)
    {
        std::cerr << "error: missing FASTA file\n";
        usage(argv[0]);
        fasta_provided = 0;
    }

    MPI_BCAST(&fasta_provided, 1, MPI_INT, root, comm);
    if (!fasta_provided) return -1;

    int fnamelen;

    if (myrank == root)
    {
        fasta_fname = argv[optind];
        fnamelen = fasta_fname.size();

        std::cout << "-DKMER_SIZE="       << KMER_SIZE       << "\n"
                  << "-DLOWER_KMER_FREQ=" << LOWER_KMER_FREQ << "\n"
                  << "-DUPPER_KMER_FREQ=" << UPPER_KMER_FREQ << "\n"
        #ifdef USE_BLOOM
                  << "-DUSE_BLOOM\n"
        #endif
                  << "\n"
                  << "int mat = "          << mat          <<   ";\n"
                  << "int mis = "          << mis          <<   ";\n"
                  << "int gap = "          << gap          <<   ";\n"
                  << "int xdrop_cutoff = " << xdrop_cutoff <<   ";\n"
                  << "String fname = \""   << fasta_fname  << "\";\n"
                  << std::endl;
    }

    MPI_BCAST(&fnamelen, 1, MPI_INT, root, comm);

    if (myrank != root) fasta_fname.assign(fnamelen, '\0');

    MPI_BCAST(fasta_fname.data(), fnamelen, MPI_CHAR, root, comm);

    return 0;
}

//  if (ret < 0)
//  {
//    parops->teardown_parallelism();
//    return ret;
//  }
//
//  /*! bcast job id */
//  // Sender
//  int job_name_length = job_name.length();
//
//  if(parops->world_proc_rank == 0)
//  {
//    MPI_Bcast(&job_name[0], job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
//  }
//  else
//  {
//    char* buf = new char[job_name_length+1];
//    buf[job_name_length] = 0;
//
//    MPI_Bcast(buf, job_name_length, MPI_CHAR, 0, MPI_COMM_WORLD);
//    std::string s(buf);
//    job_name = s;
//    delete [] buf;
//  }
//
//    int nthreads = 1;
//#ifdef THREADED
//#pragma omp parallel
//    {
//        nthreads = omp_get_num_threads();
//    }
//#endif
//
//    int nprocs, myrank;
//    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
//    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
//
//    if(myrank == 0)
//    {
//        std::cout << "Process Grid (p x p x t): " << sqrt(nprocs) << " x " << sqrt(nprocs) << " x " << nthreads << std::endl;
//    }
//
//  proc_log_file = job_name + "_rank_" + std::to_string(parops->world_proc_rank) + "_log.txt";
//  proc_log_stream.open(proc_log_file);
//
//  is_print_rank = (parops->world_proc_rank == 0);
//  std::shared_ptr<TimePod> tp = std::make_shared<TimePod>();
//  TraceUtils tu(is_print_rank);
//
//  /*! Print start time information */
//  tp->times["StartMain"] = std::chrono::system_clock::now();
//  std::time_t start_prog_time = std::chrono::system_clock::to_time_t(
//  tp->times["StartMain"]);
//
//  print_str = "\nINFO: Program started on ";
//  print_str.append(std::ctime(&start_prog_time));
//  print_str.append("\nINFO: Job ID ").append(job_name).append("\n");
//  pretty_print_config(print_str);
//  tu.print_str(print_str);

//  //////////////////////////////////////////////////////////////////////////////////////
//  // VARIALBES                                                                        //
//  //////////////////////////////////////////////////////////////////////////////////////
//
//  std::shared_ptr<DistributedFastaData> dfd;
//
//  PSpMat<PosInRead>::MPI_DCCols *Amat, *ATmat;
//  PSpMat<elba::CommonKmers>::MPI_DCCols *Bmat;
//  PSpMat<ReadOverlap>::MPI_DCCols *Rmat;
//
//  //////////////////////////////////////////////////////////////////////////////////////
//  // PIPELINE                                                                         //
//  //////////////////////////////////////////////////////////////////////////////////////
//
//  /* allocates dfd  (shared ptr) */
//  dfd = ParallelFastaParser(input_file.c_str(), idx_map_file.c_str(), parops, tp, tu);
//
//  /* allocates Amat
//   * allocates ATmat */
//  GenerateKmerByReadMatrix(dfd, Amat, ATmat, parops, tp, tu);
//
//  /* allocates Bmat
//   * deletes Amat
//   * deletes ATmat */
//  OverlapDetection(dfd, Bmat, Amat, ATmat, tp, tu);
//
//  /* allocates Rmat
//   * deletes Bmat */
//  PairwiseAlignment(dfd, Bmat, Rmat, parops, tp, tu);
//
//  //////////////////////////////////////////////////////////////////////////////////////
//  // TRANSITIVE REDUCTION                                                             //
//  //////////////////////////////////////////////////////////////////////////////////////
//
//  tp->times["StartMain:TransitiveReduction()"] = std::chrono::system_clock::now();
//
//  bool transitive_reduction = true; // use in development only
//  if (transitive_reduction)
//  {
//    TransitiveReduction(*Rmat, tu);
//  }
//
//  Rmat->Apply(Tupleize());
//
//  tp->times["EndMain:TransitiveReduction()"] = std::chrono::system_clock::now();
//
//  if(is_print_rank)
//  {
//    std::cout << "Transitive Reduction is done and okay!" << std::endl;
//  }
//
//  ////////////////////////////////////////////////////////////////////////////////////////
//  //// CONTIG EXTRACTION                                                                //
//  ////////////////////////////////////////////////////////////////////////////////////////
//
//  tp->times["StartMain:ExtractContig()"] = std::chrono::system_clock::now();
//
//  std::vector<std::string> myContigSet;
//  bool contigging = true;
//
//  if(contigging)
//  {
//    myContigSet = CreateContig(*Rmat, dfd, myoutput, tp, tu);
//  }
//
//  tp->times["EndMain:ExtractContig()"] = std::chrono::system_clock::now();
//
//  delete Rmat;
//
//  tp->times["StartMain:WriteContigs()"] = std::chrono::system_clock::now();
//
//  int64_t number_of_contigs = myContigSet.size();
//  int64_t contigs_offset = 0;
//  MPI_Exscan(&number_of_contigs, &contigs_offset, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
//
//  std::stringstream contig_filecontents;
//
//  for (int i = 0; i < myContigSet.size(); ++i)
//    contig_filecontents << ">contig" << i+contigs_offset << "\tmyrank=" << myrank << "\tmyoffset=" << i << "\n" << myContigSet[i] << "\n";
//
//  MPI_File cfh;
//  MPI_File_open(MPI_COMM_WORLD, "elba.contigs.fa", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &cfh);
//
//  std::string cfs = contig_filecontents.str();
//  const char *strout = cfs.c_str();
//
//  MPI_Offset count = strlen(strout);
//  MPI_File_write_ordered(cfh, strout, count, MPI_CHAR, MPI_STATUS_IGNORE);
//  MPI_File_close(&cfh);
//
//  MPI_Barrier(MPI_COMM_WORLD);
//  tp->times["EndMain:WriteContigs()"] = std::chrono::system_clock::now();
//
//  tp->times["EndMain"] = std::chrono::system_clock::now();
//
//  std::time_t end_prog_time = std::chrono::system_clock::to_time_t(tp->times["EndMain"]);
//
//  print_str = "INFO: Program ended on ";
//  print_str.append(std::ctime(&end_prog_time));
//  tu.print_str(print_str);
//  tu.print_str(tp->to_string());
//
//  proc_log_stream.close();
//
//    MPI_Barrier(MPI_COMM_WORLD);
//
//    return 0;
//}

//////////////////////////////////////////////////////////////////////////////////////
// INPUT COMMAND LINE PARSER                                                        //
//////////////////////////////////////////////////////////////////////////////////////

//int parse_args(int argc, char **argv)
//{
//  cxxopts::Options options("diBELLA",
//                           "Distributed Long Read to Long Read Alignment");
//
//  options.add_options()
//    (CMD_OPTION_INPUT, CMD_OPTION_DESCRIPTION_INPUT,
//     cxxopts::value<std::string>())
//    (CMD_OPTION_INPUT_SEQ_COUNT, CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT,
//     cxxopts::value<int>())
//    (CMD_OPTION_INPUT_OVERLAP, CMD_OPTION_DESCRIPTION_INPUT_OVERLAP,
//     cxxopts::value<uint64_t>())
//    (CMD_OPTION_SEED_COUNT, CMD_OPTION_DESCRIPTION_SEED_COUNT,
//     cxxopts::value<int>())
//    (CMD_OPTION_MATCH, CMD_OPTION_DESCRIPTION_MATCH,
//    cxxopts::value<int>())
//    (CMD_OPTION_MISMATCH, CMD_OPTION_DESCRIPTION_MISMATCH,
//    cxxopts::value<int>())
//    (CMD_OPTION_GAP_OPEN, CMD_OPTION_DESCRIPTION_GAP_OPEN,
//     cxxopts::value<int>())
//    (CMD_OPTION_GAP_EXT, CMD_OPTION_DESCRIPTION_GAP_EXT,
//     cxxopts::value<int>())
//    (CMD_OPTION_KMER_LENGTH, CMD_OPTION_DESCRIPTION_KMER_LENGTH,
//     cxxopts::value<int>())
//    (CMD_OPTION_KMER_STRIDE, CMD_OPTION_DESCRIPTION_KMER_STRID,
//     cxxopts::value<int>())
//    (CMD_OPTION_OVERLAP_FILE, CMD_OPTION_DESCRIPTION_OVERLAP_FILE,
//     cxxopts::value<std::string>())
//    (CMD_OPTION_ALIGN_FILE, CMD_OPTION_DESCRIPTION_ALIGN_FILE,
//     cxxopts::value<std::string>())
//    (CMD_OPTION_NO_ALIGN, CMD_OPTION_DESCRIPTION_NO_ALIGN)
//    (CMD_OPTION_FULL_ALIGN, CMD_OPTION_DESCRIPTION_FULL_ALIGN)
//    (CMD_OPTION_XDROP_ALIGN, CMD_OPTION_DESCRIPTION_XDROP_ALIGN,
//     cxxopts::value<int>())
//    (CMD_OPTION_IDX_MAP, CMD_OPTION_DESCRIPTION_IDX_MAP,
//     cxxopts::value<std::string>())
//    (CMD_OPTION_ALPH, CMD_OPTION_DESCRIPTION_ALPH,
//     cxxopts::value<std::string>())
//    (CMD_OPTION_JOB_NAME_PREFIX, CMD_OPTION_DESCRIPTION_JOB_NAME_PREFIX,
//     cxxopts::value<std::string>())
//    (CMD_OPTION_SUBS, CMD_OPTION_DESCRIPTION_SUBS,
//     cxxopts::value<int>())
//    (CMD_OPTION_LOG_FREQ, CMD_OPTION_DESCRIPTION_LOG_FREQ,
//     cxxopts::value<int>())
//    (CMD_OPTION_AF_FREQ, CMD_OPTION_DESCRIPTION_AF_FREQ,
//     cxxopts::value<int>());
//
//  auto result = options.parse(argc, argv);
//
//  bool is_world_rank0 = parops->world_proc_rank == 0;
//  if (result.count(CMD_OPTION_INPUT)) {
//    input_file = result[CMD_OPTION_INPUT].as<std::string>();
//  } else {
//    if (is_world_rank0) {
//      std::cout << "ERROR: Input file not specified" << std::endl;
//    }
//    return -1;
//  }
//
//    if (result.count(CMD_OPTION_IDX_MAP)) {
//        idx_map_file = result[CMD_OPTION_IDX_MAP].as<std::string>();
//    } else {
//        if (is_world_rank0) {
//            std::cout << "ERROR: Index map file not specified" << std::endl;
//        }
//        return -1;
//    }
//
//  // TODO - fix option parsing
//  if (result.count(CMD_OPTION_INPUT_SEQ_COUNT)) {
//    seq_count = result[CMD_OPTION_INPUT_SEQ_COUNT].as<int>();
//  } else {
//    if (is_world_rank0) {
//      std::cout << "ERROR: Input sequence count not specified" << std::endl;
//    }
//    return -1;
//  }
//
//  if (result.count(CMD_OPTION_INPUT_OVERLAP)) {
//    input_overlap = result[CMD_OPTION_INPUT_OVERLAP].as<uint64_t>();
//  } else {
//    input_overlap = 10000;
//  }
//
//  if (result.count(CMD_OPTION_SEED_COUNT)) {
//    seed_count = result[CMD_OPTION_SEED_COUNT].as<int>();
//  } else {
//    seed_count = 2;
//  }
//
//  if (result.count(CMD_OPTION_MATCH)) {
//    match = result[CMD_OPTION_MATCH].as<int>();
//  } else {
//    match = 1;
//  }
//
//  if (result.count(CMD_OPTION_MISMATCH)) {
//     mismatch_sc = result[CMD_OPTION_MISMATCH].as<int>();
//  } else {
//     mismatch_sc = -1;
//  }
//
//  if (result.count(CMD_OPTION_GAP_OPEN)) {
//    gap_open = result[CMD_OPTION_GAP_OPEN].as<int>();
//  } else {
//    gap_open = 0;
//  }
//
//  if (result.count(CMD_OPTION_GAP_EXT)) {
//    gap_ext = result[CMD_OPTION_GAP_EXT].as<int>();
//  } else {
//    gap_ext = -1;
//  }
//
//  if (result.count(CMD_OPTION_KMER_LENGTH)) {
//    klength = result[CMD_OPTION_KMER_LENGTH].as<int>();
//  } else {
//    if (is_world_rank0) {
//      std::cout << "ERROR: Kmer length not specified" << std::endl;
//    }
//    return -1;
//  }
//
//  if (result.count(CMD_OPTION_KMER_STRIDE)) {
//    kstride = result[CMD_OPTION_KMER_STRIDE].as<int>();
//  } else {
//    if (is_world_rank0) {
//      kstride = 1;
//    }
//  }
//
//  if (result.count(CMD_OPTION_OVERLAP_FILE)) {
//    write_overlaps = true;
//    overlap_file = result[CMD_OPTION_OVERLAP_FILE].as<std::string>();
//  }
//
//  if (result.count(CMD_OPTION_ALIGN_FILE)) {
//    myoutput = result[CMD_OPTION_ALIGN_FILE].as<std::string>();
//  }
//
//  if (result.count(CMD_OPTION_NO_ALIGN)) {
//    /* GGGG: You need to call the outer function anyway to calculate overhangs but the xdrop
//    case value doesn't matter, because the actual pw alignment function is not executed */
//    xdropAlign = true;
//    noAlign    = true;
//  }
//
//  if (result.count(CMD_OPTION_FULL_ALIGN)) {
//    fullAlign = true;
//    noAlign  = false;
//  }
//
//  if (result.count(CMD_OPTION_XDROP_ALIGN)) {
//    xdropAlign = true;
//    noAlign   = false;
//    xdrop = result[CMD_OPTION_XDROP_ALIGN].as<int>();
//  }
//
//  if (result.count(CMD_OPTION_ALPH)) {
//    std::string tmp = result[CMD_OPTION_ALPH].as<std::string>();
//    if (tmp == "dna"){
//      alph_t = Alphabet::DNA;
//    } else {
//      if (is_world_rank0) {
//        std::cout << "ERROR: Unsupported alphabet type " << tmp << std::endl;
//      }
//    }
//  } else {
//    alph_t = Alphabet::DNA;
//  }
//
//  if (result.count(CMD_OPTION_LOG_FREQ)) {
//    log_freq = result[CMD_OPTION_LOG_FREQ].as<int>();
//  }
//
//  if (result.count(CMD_OPTION_AF_FREQ)) {
//    afreq = result[CMD_OPTION_AF_FREQ].as<int>();
//  }
//
//  return 0;
//}
//
//auto bool_to_str = [](bool b){
//  return b ? "True" : "False";
//};
//
///*! GGGG: add input match/ mismatch_sc */
//void pretty_print_config(std::string &append_to) {
//  std::vector<std::string> params = {
//    "Input file (-i)",
//    "Original sequence count (-c)",
//    "Kmer length (k)",
//    "Kmer stride (s)",
//    "Overlap in bytes (-O)",
//    "Max seed count (--sc)",
//    "Base match score (--ma)",
//    "Base mismatch score (--mi)",
//    "Gap open penalty (-g)",
//    "Gap extension penalty (-e)",
//    "Overlap file (--of)",
//    "Alignment file (--af)",
//    "Alignment write frequency (--afreq)",
//    "No align (--na)",
//    "Full align (--fa)",
//    "Xdrop align (--xa)",
//    "Index map (--idxmap)",
//    "Alphabet (--alph)"
//  };
//
//  std::vector<std::string> vals = {
//    input_file,
//    std::to_string(seq_count),
//    std::to_string(klength),
//    std::to_string(kstride),
//    std::to_string(input_overlap),
//    std::to_string(seed_count),
//    std::to_string(match),
//    std::to_string(mismatch_sc),
//    std::to_string(gap_open),
//    std::to_string(gap_ext),
//    !overlap_file.empty() ? overlap_file  : "None",
//    !myoutput.empty()     ? myoutput    : "None",
//    !myoutput.empty()     ? std::to_string(afreq) : "None",
//    bool_to_str(noAlign),
//    bool_to_str(fullAlign),
//    bool_to_str(xdropAlign)  + (xdropAlign  ? " | xdrop: " + std::to_string(xdrop) : ""),
//    !idx_map_file.empty() ? idx_map_file : "None",
//    std::to_string(alph_t)
//  };
//
//  ushort max_length = 0;
//  for (const auto &param : params)
//  {
//    if (param.size() > max_length)
//    {
//      max_length = static_cast<ushort>(param.size());
//    }
//  }
//
//  std::string prefix = "  ";
//  append_to.append("Parameters:\n");
//
//  for (ushort i = 0; i < params.size(); ++i)
//  {
//    std::string param = params[i];
//    append_to.append(prefix).append(param).append(": ")
//      .append(get_padding(
//        static_cast<ushort>(max_length - param.size()), ""))
//      .append(vals[i]).append("\n");
//  }
//}
//
//std::string get_padding(ushort count, std::string prefix) {
//  std::string pad = std::move(prefix);
//  for (ushort i = 0; i < count; ++i) {
//    pad += " ";
//  }
//  return pad;
//}
//
//std::shared_ptr<DistributedFastaData>
//ParallelFastaParser(const char *input_file, const char *idx_map_file, const std::shared_ptr<ParallelOps>& parops, const std::shared_ptr<TimePod>& tp, TraceUtils& tu)
//{
//    tp->times["StartMain:newDFD()"] = std::chrono::system_clock::now();
//
//    std::shared_ptr<DistributedFastaData> dfd = std::make_shared<DistributedFastaData>(
//        input_file, idx_map_file, input_overlap, klength, parops, tp, tu);
//
//    tp->times["EndMain:newDFD()"] = std::chrono::system_clock::now();
//
//    if (dfd->global_count() != seq_count)
//    {
//        uint64_t final_seq_count = dfd->global_count();
//        print_str = "\nINFO: Modified sequence count\n";
//        print_str.append("  Final sequence count: ")
//                 .append(std::to_string(final_seq_count))
//                 .append(" (")
//                 .append(std::to_string((((seq_count - final_seq_count) * 100.0) / seq_count)))
//                 .append("% removed)");
//
//        seq_count = dfd->global_count();
//        print_str += "\n";
//        tu.print_str(print_str);
//    }
//    return dfd;
//}
//
//void GenerateKmerByReadMatrix(std::shared_ptr<DistributedFastaData> dfd, PSpMat<PosInRead>::MPI_DCCols*& Amat, PSpMat<PosInRead>::MPI_DCCols*& ATmat, const std::shared_ptr<ParallelOps>& parops, const std::shared_ptr<TimePod>& tp, TraceUtils& tu)
//{
//    Alphabet alph(alph_t);
//
//    tp->times["StartMain:GenerateA()"] = std::chrono::system_clock::now();
//
//    Amat = new PSpMat<PosInRead>::MPI_DCCols(elba::KmerOps::GenerateA(seq_count, dfd, klength, kstride, alph, parops, tp));
//
//    tu.print_str("Matrix A: ");
//    tu.print_str("\nLoad imbalance: " + std::to_string(Amat->LoadImbalance()) + "\n");
//
//    tp->times["EndMain:GenerateA()"] = std::chrono::system_clock::now();
//
//    Amat->PrintInfo();
//
//    ATmat = new PSpMat<PosInRead>::MPI_DCCols(*Amat);
//    tp->times["StartMain:At()"] = tp->times["EndMain:GenerateA()"];
//    ATmat->Transpose();
//    tu.print_str("Matrix At: ");
//    ATmat->PrintInfo();
//    tp->times["EndMain:At()"] = std::chrono::system_clock::now();
//}
//
//void OverlapDetection(std::shared_ptr<DistributedFastaData> dfd,
//                      PSpMat<elba::CommonKmers>::MPI_DCCols*& Bmat,
//                      PSpMat<PosInRead>::MPI_DCCols* Amat,
//                      PSpMat<PosInRead>::MPI_DCCols* ATmat,
//                      const std::shared_ptr<TimePod>& tp, TraceUtils& tu)
//{
//    tp->times["StartMain:AAt()"] = std::chrono::system_clock::now();
//
//    // @GGGG-TODO: check vector version (new one stack error)
//    Bmat = new PSpMat<elba::CommonKmers>::MPI_DCCols(Mult_AnXBn_DoubleBuff<KmerIntersectSR_t, elba::CommonKmers, PSpMat<elba::CommonKmers>::DCCols>(*Amat, *ATmat));
//
//    delete Amat;
//    delete ATmat;
//
//    tp->times["EndMain:AAt()"] = std::chrono::system_clock::now();
//
//    // @GGGG-TODO: remove proc_log_stream
//    tu.print_str(
//        "Matrix AAt: Overlaps after k-mer finding (nnz(C) - diagonal): "
//        + std::to_string(Bmat->getnnz() - seq_count)
//        + "\nLoad imbalance: " + std::to_string(Bmat->LoadImbalance()) + "\n");
//
//    tu.print_str("Matrix B, i.e AAt: ");
//    Bmat->PrintInfo();
//
//    /*! Wait until sequence exchange is complete */
//    tp->times["StartMain:DfdWait()"] = std::chrono::system_clock::now();
//    if (!dfd->is_ready())
//    {
//      dfd->wait();
//    }
//    tp->times["EndMain:DfdWait()"] = std::chrono::system_clock::now();
//}
//
//void PairwiseAlignment(std::shared_ptr<DistributedFastaData> dfd, PSpMat<elba::CommonKmers>::MPI_DCCols* Bmat, PSpMat<ReadOverlap>::MPI_DCCols*& Rmat, const std::shared_ptr<ParallelOps>& parops, const std::shared_ptr<TimePod>& tp, TraceUtils& tu)
//{
//  uint64_t n_rows, n_cols;
//  n_rows = n_cols = dfd->global_count();
//  int gr_rows = parops->grid->GetGridRows();
//  int gr_cols = parops->grid->GetGridCols();
//
//  int gr_col_idx = parops->grid->GetRankInProcRow();
//  int gr_row_idx = parops->grid->GetRankInProcCol();
//
//  uint64_t avg_rows_in_grid = n_rows / gr_rows;
//  uint64_t avg_cols_in_grid = n_cols / gr_cols;
//  uint64_t row_offset = gr_row_idx * avg_rows_in_grid;  // first row in this process
//  uint64_t col_offset = gr_col_idx * avg_cols_in_grid;	// first col in this process
//
//  DistributedPairwiseRunner dpr(dfd, Bmat->seqptr(), Bmat, afreq, row_offset, col_offset, parops);
//
//  double mytime = MPI_Wtime();
//  tp->times["StartMain:DprAlign()"] = std::chrono::system_clock::now();
//  ScoringScheme scoring_scheme(match, mismatch_sc, gap_ext);
//
//  PairwiseFunction* pf = nullptr;
//  uint64_t local_alignments = 1;
//
//  // Output intermediate matrix post-alignment
//  //std::string candidatem = myoutput;
//  //candidatem += ".candidatematrix.mm";
//  //Bmat->ParallelWriteMM(candidatem, true, elba::CkOutputMMHandler());
//
//  if(xdropAlign)
//  {
//    pf = new SeedExtendXdrop (scoring_scheme, klength, xdrop, seed_count);
//    dpr.run_batch(pf, proc_log_stream, log_freq, ckthr, aln_score_thr, tu, noAlign, klength, seq_count);
//	  local_alignments = static_cast<SeedExtendXdrop*>(pf)->nalignments;
//  }
//  else if(fullAlign)
//  {
//    pf = new FullAligner(scoring_scheme);
//    dpr.run_batch(pf, proc_log_stream, log_freq, ckthr, aln_score_thr, tu, noAlign, klength, seq_count);
//	  local_alignments = static_cast<FullAligner*>(pf)->nalignments;
//  }
//
//  tp->times["EndMain:DprAlign()"] = std::chrono::system_clock::now();
//  delete pf;
//
//  uint64_t total_alignments = 0;
//  MPI_Reduce(&local_alignments, &total_alignments, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
//
//  // total_alignments should be zero if "noAlign" is true
//  if(is_print_rank)
//  {
//    std::cout << "Final alignment (L+U-D) count: " << 2 * total_alignments << std::endl;
//  }
//
//  // Output intermediate matrix post-alignment
//  //std::string postalignment = myoutput;
//  //postalignment += ".resultmatrix.mm";
//  //Bmat->ParallelWriteMM(postalignment, true, elba::CkOutputHandler());
//
//  Rmat = new PSpMat<ReadOverlap>::MPI_DCCols(*Bmat);
//
//  //Rmat->ParallelWriteMM("alignment.mtx", true, ReadOverlapGraphHandler());
//
//  delete Bmat;
//}
