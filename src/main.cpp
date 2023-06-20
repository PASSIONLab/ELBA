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
#include "TransitiveReduction.hpp"
#include "ContigGeneration.hpp"
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
CT<int64_t>::PDistVec find_contained_reads(const CT<Overlap>::PSpParMat& R);
std::string get_overlap_paf_name();
std::string get_string_paf_name();

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
        std::unique_ptr<CT<Overlap>::PSpParMat> R, S;
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

        /*
         * The next steps can be understood by first describing what @kmermap
         * is. @kmermap is an unordered_map (basically an associative container)
         * mapping k-mers to "k-mer count entries". Formally, a "k-mer count entry"
         * is a triple (READIDS, POSITIONS, count) where READIDS is an array of global
         * read IDs, POSITIONS is an array of read positions, and count is the number
         * of distinct times that k-mer has been found in the FASTA sequences (globally).
         *
         * It should be noted that READIDS and POSITIONS are arrays with a fixed size,
         * determined by the compile-time parameter UPPER_KMER_FREQ. This is because
         * we only want to store k-mers that appear <= UPPER_KMER_FREQ different times
         * in the input.
         *
         * @kmermap is a "distributed" hash table in the sense that each processor
         * has its own local instance which is individually responsible for a different
         * partition of k-mers. That is to say: given a k-mer @s, there is exactly
         * one processor that is responsible for storing @s and its associated
         * k-mer count entry. This processor is determined by hashing @s.
         *
         * So what does get_kmer_count_map_keys() do exactly? Basically it does
         * 4 things:
         *
         *    1. Globally estimates the number of distinct k-mers in the dataset using
         *       the HyperLogLog data structure. This estimate is used to allocate
         *       memory on each process for its local partition of the distributed hash table.
         *
         *    2. Each process in parallel computes all the k-mers (not required to be distinct)
         *       that exist in its local FASTA partition (@mydna, which was obtained by index.getmydna()).
         *
         *    3. All processors collectively send their locally found k-mers to their proper
         *       destinations. In parallel, each processor receives incoming k-mers assigned to
         *       it, and filters out likely singletons using a Bloom filter approach. Only
         *       distinct k-mers are kept.
         *
         *    4. The received distinct k-mers are used to initialize @kmermap with empty "kmer count
         *       entries." Those entries are filled out in the second pass implemented in
         *       @get_kmer_count_map_values().
         *
         */
        timer.start();
        kmermap = get_kmer_count_map_keys(mydna, commgrid);
        timer.stop_and_log("collecting distinct k-mers");

        /*
         * Now that every process has its local partition of the distributed k-mer hash table
         * initialized with all the keys (k-mers) it needs, we do a second pass over every k-mer
         * and send them to their destinations, this time including information about which read ID
         * the k-mer was parsed from, and the position of the k-mer within that read. Using the
         * Bloom filter, we can quickly query received k-mers and discard those k-mers which we know
         * aren't keys in the local @kmermap.
         *
         * For each received k-mer that passes through the Bloom filter (is accepted), its
         * corresponding triple (READIDS, POSITIONS, count) is updated. This way, the local
         * process is able to quickly find all the reads (via their global IDs) that contain a particular
         * k-mer, and quickly find the position within that read where the k-mer is located. The
         * count parameter merely states how many times that k-mer has been found in the dataset,
         * and is therefore equivalent to the number of used entries of READIDS and POSITIONS.
         *
         * Suppose that a k-mer @s appears more that UPPER_KMER_FREQ times in the input. Then
         * it is guaranteed that, eventually, the processor responsible for storing @s will
         * receive an instance of @s (plus a global read id and position where @s came from) that
         * will put it over the capacity of POSITIONS and READIDS. We therefore always check
         * if a received k-mer will push us over this threshold, and if it does, we DELETE the
         * k-mer key of @s on the owner processor. That way, any other instances of @s from
         * other reads are discarded because we only record entries for k-mers that exist in
         * the hash table.
         *
         * Finally, once the collective exchange is finished, we delete all k-mer keys of
         * k-mers that appeared less than LOWER_KMER_FREQ times. The result is a distributed
         * hash table mapping reliable k-mers (k-mers that appear within the defined frequency bounds)
         * to their corresponding k-mer count entries.
         */
        timer.start();
        get_kmer_count_map_values(mydna, *kmermap, commgrid);
        timer.stop_and_log("counting recording k-mer seeds");

        print_kmer_histogram(*kmermap, commgrid);

        /*
         * Now that all the reliable k-mers and their locations have been computed and stored
         * in the distributed hash table @kmermap, it is time to construct the distributed k-mer sparse
         * matrix @A. The purpose of @A is to facilitate read overlap detection via an SpGEMM
         * operation using a custom semiring. More on that will be discussed later. For now,
         * this is what @A is:
         *
         *    Let M = number of reads in FASTA;
         *    Let N = number of distinct k-mers (keys) currently stored in distributed k-mer hash table;
         *    Let L = total number of k-mer instances (values) currently stored in the distributed k-mer hash table;
         *
         *    Then @A is an M-by-N distributed sparse matrix with L nonzeros, where a nonzero at
         *    @A(i,j) represents an instance of a k-mer (with global id j) found in
         *    read sequence i (global id). The global k-mer ids are computed using a prefix
         *    scan of the stored k-mer keys in the distributed hash table. The actual
         *    value stored by the nonzero is the POSITION of k-mer j within read i.
         *
         * A few quick observations on what this means:
         *
         *    The number of nonzeros in the row @A(i,:) is the number of distinct reliable k-mers
         *    found within the sequence with global id i.
         *
         *    The number of nonzeros in the column @A(:,j) is the number of different sequences
         *    that contain the reliable k-mer with id j as a subsequence.
         *
         * Other similar observations about the nature of @A can be made, but hopefully it
         * is clear by now what @A is.
         */
        timer.start();
        A = create_kmer_matrix(mydna, *kmermap, commgrid);
        timer.stop_and_log("creating k-mer matrix");

        /*
         * Once @A has been constructed, we have no more use for the distributed k-mer hash table
         * so we release all its memory.
         */
        kmermap.reset();

        /*
         * The SpGEMM overlap detection phase requires both @A and its transpose @AT.
         */
        timer.start();
        AT = std::make_unique<CT<PosInRead>::PSpParMat>(*A);
        AT->Transpose();
        timer.stop_and_log("copying and transposing k-mer matrix");
        elbalog.log_kmer_matrix(*A);

        /*
         * TODO: comment this.
         */
        timer.start();
        B = create_seed_matrix(*A, *AT);
        timer.stop_and_log("creating seed matrix (spgemm)");

        A.reset();
        AT.reset();

        elbalog.log_seed_matrix(*B);

        dfd.wait();

        /*
         * In order to obtain reliable overlaps, we need to do some alignments.
         * The @B matrix provides the seeds (common k-mers) from which we
         * anchor and extend our alignments using the X-drop algorithm.
         * In an embarassingly parallel manner, we run seed-and-extend
         * alignments on each nonzero (actually only half since @B is symmetric)
         * and then prune the alignments that appear spurious.
         */
        timer.start();
        R = PairwiseAlignment(dfd, *B, mat, mis, gap, xdrop_cutoff);
        timer.stop_and_log("pairwise alignment");

        elbalog.log_overlap_matrix(*R);
        parallel_write_paf(*R, dfd, get_overlap_paf_name().c_str());

        timer.start();
        auto contained = find_contained_reads(*R);
        R->PruneFull(contained, contained);

        S = TransitiveReduction(*R);
        timer.stop_and_log("contained read removal and transitive reduction");
        R.reset();

        parallel_write_paf(*S, dfd, get_string_paf_name().c_str());
        S->ParallelWriteMM("S.mtx", true, Overlap::IOHandler());

        // std::vector<std::string> contigs = GenerateContigs(*S, mydna, dfd);

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

    for (int64_t i = 0; i < dcsc->nzc; ++i)
        for (int64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
        {
            int64_t localrow = dcsc->ir[j];
            int64_t localcol = dcsc->jc[i];
            int64_t globalrow = localrow + dfd.getrowstartid();
            int64_t globalcol = localcol + dfd.getcolstartid();

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

CT<int64_t>::PDistVec find_contained_reads(const CT<Overlap>::PSpParMat& R)
{
    CT<int>::PSpParMat containedQmat = const_cast<CT<Overlap>::PSpParMat&>(R).Prune([](const Overlap& o) { return !o.containedQ; }, false);
    CT<int>::PSpParMat containedTmat = const_cast<CT<Overlap>::PSpParMat&>(R).Prune([](const Overlap& o) { return !o.containedT; }, false);
    CT<int>::PDistVec containedQvec = containedQmat.Reduce(Row, std::logical_or<int>(), 0);
    CT<int>::PDistVec containedTvec = containedTmat.Reduce(Column, std::logical_or<int>(), 0);
    CT<int>::PDistVec contained = containedQvec;
    contained.EWiseOut(containedTvec, std::logical_or<int>(), contained);
    return contained.FindInds([](const int& v) { return v > 0; });
}

std::string get_string_paf_name()
{
    std::ostringstream ss;
    ss << output_prefix << ".string.paf";
    return ss.str();
}

std::string get_overlap_paf_name()
{
    std::ostringstream ss;
    ss << output_prefix << ".overlap.paf";
    return ss.str();
}
