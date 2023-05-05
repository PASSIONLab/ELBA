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
#include "DistributedFastaData.hpp"
#include "Kmer.hpp"
#include "KmerOps.hpp"
#include "SharedSeeds.hpp"
#include "PairwiseAlignment.hpp"
#include "MPITimer.hpp"

int returncode;
std::string fasta_fname;
std::string output_prefix = "elba";

std::string getmatfname(const std::string matname);
std::string getpafname();

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

int parse_cli(int argc, char *argv[]);
void print_kmer_histogram(const KmerCountMap& kmermap, std::shared_ptr<CommGrid> commgrid);
void parallel_write_paf(const CT<Overlap>::PSpParMat& R, DistributedFastaData& dfd, char const *pafname);

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    {
        /*
         * CombBLAS 2D processor grid uses MPI_COMM_WORLD
         * as its world communicator (comm).
         */
        std::shared_ptr<CommGrid> commgrid;

        comm = MPI_COMM_WORLD;
        commgrid.reset(new CommGrid(comm, 0, 0));
        myrank = commgrid->GetRank();
        nprocs = commgrid->GetSize();

        /*
         * Parse command-line from root process and broadcast.
         */
        if (parse_cli(argc, argv) != 0)
            goto err;

        /***********************************************************/
        /************************ START ****************************/
        /***********************************************************/

        MPITimer timer(comm);

        timer.start();
        FastaIndex index(fasta_fname, commgrid);
        timer.stop_and_log("indexing");

        timer.start();
        std::shared_ptr<DnaBuffer> mydna = index.getmydna();
        timer.stop_and_log("parsing and compressing");

        DistributedFastaData dfd(index);

        dfd.collect_sequences(mydna);

        KmerCountMap kmermap = GetKmerCountMapKeys(*mydna, commgrid);

        GetKmerCountMapValues(*mydna, kmermap, commgrid);

        print_kmer_histogram(kmermap, commgrid);

        std::unique_ptr<CT<PosInRead>::PSpParMat> A, AT;

        A = create_kmer_matrix(*mydna, kmermap, commgrid);
        AT = std::make_unique<CT<PosInRead>::PSpParMat>(*A);
        AT->Transpose();

        size_t numreads = A->getnrow();
        size_t numkmers = A->getncol();
        size_t numseeds = A->getnnz();

        if (!myrank)
        {
           std::cout << "K-mer matrix A has " << numreads << " rows (readids), " << numkmers << " columns (k-mers), and " << numseeds << " nonzeros (k-mer seeds)\n" << std::endl;
        }
        MPI_Barrier(comm);


        #if LOG_LEVEL >= 2
        A->ParallelWriteMM(getmatfname("A.mtx").c_str(), true);
        #endif


        CT<SharedSeeds>::PSpParMat B = Mult_AnXBn_DoubleBuff<SharedSeeds::Semiring, SharedSeeds, CT<SharedSeeds>::PSpDCCols>(*A, *AT);
        A.reset(), AT.reset();

        size_t numsharedseeds_before, numsharedseeds_after;

        numsharedseeds_before = B.getnnz();
        B.Prune([](const SharedSeeds& nz) { return nz.getnumshared() <= 1; });
        numsharedseeds_after = B.getnnz();

        #if LOG_LEVEL >= 1
        if (!myrank)
        {
            std::cout << "Pruned " << (numsharedseeds_before - numsharedseeds_after) << " nonzeros (overlap seeds) with not enough support\n\n";
            std::cout << "Overlap matrix B has " << numreads << " rows (readids), " << numreads << " columns (readids), and " << numsharedseeds_after << " nonzeros (overlap seeds)\n" << std::endl;
        }
        MPI_Barrier(comm);
        #endif

        #if LOG_LEVEL >= 2
        B.ParallelWriteMM(getmatfname("B.mtx").c_str(), true, SharedSeeds::IOHandler());
        #endif

        dfd.wait();

        std::vector<std::string> names = index.bcastnames();

        CT<Overlap>::PSpParMat R = PairwiseAlignment(dfd, B, mat, mis, gap, xdrop_cutoff);

        size_t numoverlaps_before = R.getnnz();
        R.Prune([](const Overlap& nz) { return !nz.passed; });
        size_t numoverlaps_after = R.getnnz();

        #if LOG_LEVEL >= 1
        if (!myrank)
        {
            std::cout << "Pruned " << (numoverlaps_before - numoverlaps_after) << " nonzeros (overlaps) with poor alignments\n\n";
        }
        #endif

        #if LOG_LEVEL >= 2
        R.ParallelWriteMM(getmatfname("R.mtx").c_str(), true, Overlap::IOHandler());
        #endif

        parallel_write_paf(R, dfd, getpafname().c_str());

        /***********************************************************/
        /************************* END *****************************/
        /***********************************************************/

        goto done;
    }

err: returncode = -1;
done: MPI_Finalize();
      return returncode;
}

void usage(char const *prg)
{
    std::cerr << "Usage: " << prg << " [options] <reads.fa>\n"
              << "Options: -x INT   x-drop alignment threshold [" <<  xdrop_cutoff               << "]\n"
              << "         -A INT   matching score ["             <<  mat                        << "]\n"
              << "         -B INT   mismatch penalty ["           << -mis                        << "]\n"
              << "         -G INT   gap penalty ["                << -gap                        << "]\n"
              << "         -o STR   output file name prefix "     <<  std::quoted(output_prefix) << "\n"
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

        while ((c = getopt(argc, argv, "x:A:B:G:o:h")) >= 0)
        {
            if      (c == 'A') params[0] =  atoi(optarg);
            else if (c == 'B') params[1] = -atoi(optarg);
            else if (c == 'G') params[2] = -atoi(optarg);
            else if (c == 'x') params[3] =  atoi(optarg);
            else if (c == 'o') output_prefix = std::string(optarg);
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

    int fnamelen, onamelen;

    if (myrank == root)
    {
        fasta_fname = argv[optind];
        fnamelen = fasta_fname.size();
        onamelen = output_prefix.size();

        std::cout << "-DKMER_SIZE="            << KMER_SIZE            << "\n"
                  << "-DLOWER_KMER_FREQ="      << LOWER_KMER_FREQ      << "\n"
                  << "-DUPPER_KMER_FREQ="      << UPPER_KMER_FREQ      << "\n"
                  << "-DMAX_SEEDS="            << MAX_SEEDS            << "\n"
                  << "-DMPI_HAS_LARGE_COUNTS=" << MPI_HAS_LARGE_COUNTS << "\n"
        #ifdef USE_BLOOM
                  << "-DUSE_BLOOM\n"
        #endif
                  << "\n"
                  << "int mat = "              << mat                        << ";\n"
                  << "int mis = "              << mis                        << ";\n"
                  << "int gap = "              << gap                        << ";\n"
                  << "int xdrop_cutoff = "     << xdrop_cutoff               << ";\n"
                  << "String fname = "         << std::quoted(fasta_fname)   << ";\n"
                  << "String output_prefix = " << std::quoted(output_prefix) << ";\n"
                  << std::endl;
    }

    MPI_BCAST(&fnamelen, 1, MPI_INT, root, comm);
    MPI_BCAST(&onamelen, 1, MPI_INT, root, comm);

    if (myrank != root)
    {
        fasta_fname.assign(fnamelen, '\0');
        output_prefix.assign(onamelen, '\0');
    }

    MPI_BCAST(fasta_fname.data(), fnamelen, MPI_CHAR, root, comm);
    MPI_BCAST(output_prefix.data(), onamelen, MPI_CHAR, root, comm);

    return 0;
}

void print_kmer_histogram(const KmerCountMap& kmermap, std::shared_ptr<CommGrid> commgrid)
{
    int maxcount = std::accumulate(kmermap.cbegin(), kmermap.cend(), 0, [](int cur, const auto& entry) { return std::max(cur, std::get<2>(entry.second)); });

    MPI_Allreduce(MPI_IN_PLACE, &maxcount, 1, MPI_INT, MPI_MAX, commgrid->GetWorld());

    std::vector<int> histo(maxcount+1, 0);

    for (auto itr = kmermap.cbegin(); itr != kmermap.cend(); ++itr)
    {
        int cnt = std::get<2>(itr->second);
        assert(cnt >= 1);
        histo[cnt]++;
    }

    MPI_Allreduce(MPI_IN_PLACE, histo.data(), maxcount+1, MPI_INT, MPI_SUM, commgrid->GetWorld());

    int myrank = commgrid->GetRank();

    if (!myrank)
    {
        std::cout << "#count\tnumkmers" << std::endl;

        for (int i = 1; i < histo.size(); ++i)
        {
            if (histo[i] > 0)
            {
                std::cout << i << "\t" << histo[i] << std::endl;
            }
        }
        std::cout << std::endl;
    }

    MPI_Barrier(commgrid->GetWorld());
}

void parallel_write_paf(const CT<Overlap>::PSpParMat& R, DistributedFastaData& dfd, char const *pafname)
{
    auto index = dfd.getindex();
    auto commgrid = index.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();
    MPI_Comm comm = commgrid->GetWorld();

    std::vector<std::string> names = index.bcastnames();

    auto dcsc = R.seqptr()->GetDCSC();

    std::ostringstream ss;

    for (uint64_t i = 0; i < dcsc->nzc; ++i)
        for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
        {
            uint64_t localrow = dcsc->ir[j];
            uint64_t localcol = dcsc->jc[i];
            uint64_t globalrow = localrow + dfd.getrowstartid();
            uint64_t globalcol = localcol + dfd.getcolstartid();

            const Overlap& o = dcsc->numx[j];

            int maplen = std::max(std::get<0>(o.end) - std::get<0>(o.beg), std::get<1>(o.end) - std::get<1>(o.end));

            ss << names[globalrow] << "\t" << std::get<0>(o.len) << "\t" << std::get<0>(o.beg) << "\t" << std::get<0>(o.end) << "\t" << "+-"[static_cast<int>(o.rc)] << "\t"
               << names[globalcol] << "\t" << std::get<1>(o.len) << "\t" << std::get<1>(o.beg) << "\t" << std::get<1>(o.end) << "\t" << o.score << "\t" << maplen << "\n";
        }

    std::string pafcontents = ss.str();

    MPI_Count count = pafcontents.size();
    MPI_File fh;
    MPI_File_open(comm, pafname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_ordered(fh, pafcontents.c_str(), count, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}

std::string getmatfname(const std::string matname)
{
    std::ostringstream ss;
    ss << output_prefix << "." << matname;
    return ss.str();
}

std::string getpafname()
{
    std::ostringstream ss;
    ss << output_prefix << ".paf";
    return ss.str();
}
