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
#include <iterator>
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
#include "PruneChimeras.hpp"
#include "MPITimer.hpp"
#include "ELBALogger.hpp"

int returncode;
std::string fasta_fname;
std::string output_prefix = "elba";

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
std::string getpafname();

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

        std::unique_ptr<CT<PosInRead>::PSpParMat> A, AT;
        std::unique_ptr<CT<SharedSeeds>::PSpParMat> B;
        std::unique_ptr<CT<Overlap>::PSpParMat> R;
        std::unique_ptr<KmerCountMap> kmermap;

        std::ostringstream ss;
        ELBALogger elbalog(output_prefix, comm);
        MPITimer timer(comm), walltimer(comm);
        walltimer.start();

        /*
         * FastaIndex @index is the structure responsible for reading
         * the .fai index file and telling each process which read sequences
         * it is responsible for parsing, compressing and storing.
         */
        timer.start();
        FastaIndex index(fasta_fname, commgrid);
        ss << "reading " << std::quoted(index.get_faidx_fname()) << " and scattering to all MPI tasks";
        timer.stop_and_log(ss.str().c_str());
        ss.clear(); ss.str("");

        /*
         * DnaBuffer @mydna stores the compressed read sequences assigned
         * to it by the .fai index file, as determined by @index.
         */
        timer.start();
        DnaBuffer mydna = index.getmydna();
        ss << "reading and 2-bit encoding " << std::quoted(index.get_fasta_fname()) << " sequences in parallel";
        timer.stop_and_log(ss.str().c_str());
        ss.clear(); ss.str("");

        /*
         * DistributedFastaData @dfd is the structure responsible for telling
         * each process which read sequences it will need according to the
         * 2D processor grid. The constructor determines this from the FastaIndex
         * @index.
         */
        DistributedFastaData dfd(index);

        /*
         * Initiate non-blocking send and receive calls. The end goal
         * is that every processor has the reads it requires according
         * the original construction above. These are stored as DnaBuffer
         * objects owned by @dfd.
         *
         * Because this is non-blocking, we have to call dfd.wait() to
         * make sure every processor has finished its sends and receives.
         * We do this after the k-mer counting step, because that step
         * only needs access to the sequences stored in mydna.
         */
        dfd.collect_sequences(mydna);

        timer.start();
        kmermap = get_kmer_count_map_keys(mydna, commgrid);
        timer.stop_and_log("collecting distinct k-mers");

        timer.start();
        get_kmer_count_map_values(mydna, *kmermap, commgrid);
        timer.stop_and_log("counting recording k-mer seeds");

        print_kmer_histogram(*kmermap, commgrid);

        timer.start();
        A = create_kmer_matrix(mydna, *kmermap, commgrid);
        timer.stop_and_log("creating k-mer matrix");

        kmermap.reset();

        timer.start();
        AT = std::make_unique<CT<PosInRead>::PSpParMat>(*A);
        AT->Transpose();
        timer.stop_and_log("copying and transposing k-mer matrix");

        elbalog.log_kmer_matrix(*A);

        timer.start();
        B = create_seed_matrix(*A, *AT);
        timer.stop_and_log("creating seed matrix (spgemm)");

        A.reset();
        AT.reset();

        elbalog.log_seed_matrix(*B);

        dfd.wait();

        timer.start();
        R = PairwiseAlignment(dfd, *B, mat, mis, gap, xdrop_cutoff);
        timer.stop_and_log("pairwise alignment");

        elbalog.log_overlap_matrix(*R);

        parallel_write_paf(*R, dfd, getpafname().c_str());

        R.reset();
        walltimer.stop_and_log("wallclock");

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
                  << "String output_prefix = " << std::quoted(output_prefix) << ";\n\n"
                  << "MPI processes = " << nprocs << "\n"
                  << "rows/columns in 2D processor grid = " << static_cast<int>(std::sqrt(nprocs)) << "\n"
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
    #if LOG_LEVEL >= 2
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
    #endif
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
               << names[globalcol] << "\t" << std::get<1>(o.len) << "\t" << std::get<1>(o.beg) << "\t" << std::get<1>(o.end) << "\t" << o.score << "\t" << maplen << "\t255\n";
        }

    std::string pafcontents = ss.str();

    MPI_Count count = pafcontents.size();
    MPI_File fh;
    MPI_File_open(comm, pafname, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    MPI_File_write_ordered(fh, pafcontents.c_str(), count, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&fh);
}

std::string getpafname()
{
    std::ostringstream ss;
    ss << output_prefix << ".paf";
    return ss.str();
}
