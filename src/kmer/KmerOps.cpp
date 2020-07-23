// Created by Saliya Ekanayake on 10/15/19.

#include <numeric>
#include "../../include/kmer/KmerOps.hpp"
#include "../../include/NearestKmers2.hpp"

/*! GGGG: define this type somewhere */
KmerCountsType *kmercounts = NULL;

namespace dibella {
  PSpMat<MatrixEntry>::MPI_DCCols KmerOps::generate_A(uint64_t seq_count,
      std::shared_ptr<DistributedFastaData> &dfd, ushort k, ushort s,
      Alphabet &alph, const std::shared_ptr<ParallelOps> &parops,
      const std::shared_ptr<TimePod> &tp, std::unordered_set<Kmer, Kmer>& local_kmers) {

    char *buff;
    ushort len;
    uint64_t start_offset, end_offset_inclusive;
    uvec_64 lrow_ids, lcol_ids;
    std::vector<MatrixEntry> lvals;
    uint64_t offset = dfd->global_start_idx();
    FastaData *lfd = dfd->lfd();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: Cardinality estimation without heavy hitters
/////////////////////////////////////////////////////////////////////////////////////////////////////////

    tp->times["start_kmerop:gen_A:loop_add_kmers()"] = std::chrono::system_clock::now();

    /*! GGGG: cardinality estimate */
    double tstart = MPI_Wtime(); 
    unsigned int readsxproc = 0;

    int64_t sums[3] = {0, 0, 0};
    int64_t &cardinality = sums[0];
    int64_t &totreads    = sums[1];
    int64_t &totbases    = sums[2];

    totreads = lfd->local_count();
    for (uint64_t lseq_idx = 0; lseq_idx < lfd->local_count(); ++lseq_idx)
    {
      /*! GGGG: loading sequence string in buff */
      buff = lfd->get_sequence(lseq_idx, len, start_offset, end_offset_inclusive);
      totbases += len;
    }

    if (totreads > 0)
    {
      /*! GGGG: double check k == kmer len */
      int64_t kmersxread = ((totbases + totreads - 1) / totreads) - k + 1;
      cardinality += kmersxread * totreads;
    }

    MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);

    if (myrank == 0)
    {
      std::cout << "Estimated cardinality: " << cardinality << " totreads: " << totreads << " totbases: " << totbases << std::endl;
    }

    /*! GGGG: from dibella v1 this is baseline for 10M kmers */
    if (cardinality < 10000000) cardinality = 10000000;

    /*! Assume maximum of 90% of the kmers are unique, because at least some have to repeat */
    cardinality = 0.9 * cardinality / (double) nprocs;
    readsxproc = totreads / nprocs;

    double tcardinalitye = MPI_Wtime() - tstart;

    /*! GGGG: print using PASTIS way */
    if (myrank == 0)
    {
      std::cout << "Estimated cardinality: " << cardinality << " totreads: " << totreads << " totbases: " << totbases << std::endl;
    }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: First pass k-mer counter
/////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*! Prepare for first pass over input data with k-mer extraction */
    ASSERT(kmercounts == NULL, "");

    /*! GGGG: what is this? */
    kmercounts = new KmerCountsType();

    /*! Assume we have at least 10x depth of kmers to store and a minimum 64k */
    uint64_t reserve = cardinality * 0.1 + 64000; 
    /*! This is the maximum supported by khash */
    if (reserve >= 4294967291ul) reserve = 4294967290ul; // 

    /*! GGGG: define this macros */
#ifdef KHASH
    LOGF("Reserving %lld entries in KHASH for cardinality %lld\n", (lld) reserve, (lld) cardinality);
    kmercounts->reserve(reserve);
#else
    LOGF("Reserving %lld entries in VectorMap for cardinality %lld\n", (lld) reserve, (lld) cardinality);
    kmercounts->reserve(reserve);
#endif
    DBG("Reserved kmercounts\n");

    /* Initialize readNameMap for storing ReadID -> names/tags of reads */
    /* GGGG: define ReadId type */
    std::unordered_map<ReadId, std::string>* readNameMap = new std::unordered_map<ReadId, std::string>();

    /* First pass */
    ReadId myReadStartIndex = 0;

    /*! GGGG: let's extract the function (I'll separate later once I understood what's going on) */
    // int nReads = ProcessFiles(allfiles, 1, cardinality, cached_io, base_dir, myReadStartIndex, *readNameMap); // determine final hash-table entries using bloom filter
    
    /////////////////////////////////////////////
    // This is gonna be a stand-alone function //
    /////////////////////////////////////////////

    /*! GGGG: KmerCounter.cpp */

    /////////////////////////////////////////////
    // End of function                         //
    /////////////////////////////////////////////   
    

    DBG("my nreads=%lld\n", nReads);

    time_first_data_pass = MPI_Wtime() - time_temp;
    serial_printf("%s: 1st input data pass, elapsed time: %0.3f s\n", __FUNCTION__, time_first_data_pass);
    time_temp = MPI_Wtime();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: Second pass k-mer counter
/////////////////////////////////////////////////////////////////////////////////////////////////////////

    /*! GGGG: prepare for second pass over the input data, extracting k-mers with read ID's, names, etc. */
    /* Exchange number of reads-per-processor to calculate read indices */
    uint64_t sndReadCounts[nprocs];

    for (int i = 0; i < nprocs; i++)
    {
      sndReadCounts[i] = nReads;
    }

    uint64_t recvReadCounts[nprocs];

    CHECK_MPI(MPI_Alltoall(sndReadCounts, 1, MPI_UINT64_T, recvReadCounts, 1, MPI_UINT64_T, MPI_COMM_WORLD));
    
    myReadStartIndex = 1;
    for (int i = 0; i < myrank; i++)
    {
      myReadStartIndex += recvReadCounts[i];
    }

#ifdef DEBUG
      if(myrank == nprocs - 1)
        for (int i = 0; i < nprocs; i++)
          LOGF("recvReadCounts[%lld] = %lld\n", recvReadCounts[i]);
#endif
    DBG("my startReadIndex = %lld\n", myReadStartIndex);

    CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));

    /* Initialize end read ranges */
    ReadId readRanges[nprocs];
    ReadId last = 0;

    for (int i = 0; i < nprocs; i++)
    {
      readRanges[i] = last + recvReadCounts[i];
      last += recvReadCounts[i];
    }

    DBG("My read range is [%lld - %lld]\n", (myrank==0? 1 : readRanges[myrank-1]+1), readRanges[myrank]);

    /* Second pass */
    ProcessFiles(allfiles, 2, cardinality, cached_io, base_dir, myReadStartIndex, *readNameMap);

    time_second_data_pass = MPI_Wtime() - time_temp;
    serial_printf("%s: 2nd input data pass, elapsed time: %0.3f s\n", __FUNCTION__, time_second_data_pass);
    time_temp = MPI_Wtime();

    int64_t sendbuf = kmercounts->size(); 
    int64_t recvbuf, totcount, maxcount;

    CHECK_MPI(MPI_Exscan(&sendbuf, &recvbuf, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
    CHECK_MPI(MPI_Allreduce(&sendbuf, &totcount, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD));
    CHECK_MPI(MPI_Allreduce(&sendbuf, &maxcount, 1, MPI_LONG_LONG, MPI_MAX, MPI_COMM_WORLD));

    int64_t totkmersprocessed;

    CHECK_MPI(MPI_Reduce(&kmersprocessed, &totkmersprocessed, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD));

    time_computing_loadimbalance = MPI_Wtime() - time_temp;

    if(myrank  == 0)
    {
      cout << __FUNCTION__ << ": Total number of stored k-mers: " << totcount << endl;

      double imbalance = static_cast<double>(nprocs * maxcount) / static_cast<double>(totcount);  
      
      cout << __FUNCTION__ << ": Load imbalance for final k-mer counts: " << imbalance << endl;
      cout << __FUNCTION__ << ": CardinalityEstimate " << static_cast<double>(totkmersprocessed) / (MEGA * max((time_cardinality_est),0.001) * nprocs) << " MEGA k-mers per sec/proc in " << (time_cardinality_est) << " seconds"<< endl;
      cout << __FUNCTION__ << ": Bloom filter + hash table (key) initialization " << static_cast<double>(totkmersprocessed) / (MEGA * max((time_first_data_pass),0.001) * nprocs) << " MEGA k-mers per sec/proc in " << (time_first_data_pass) << " seconds" << endl;
      cout << __FUNCTION__ << ": Hash table (value) initialization  " << static_cast<double>(totkmersprocessed) / (MEGA * max((time_second_data_pass),0.001) * nprocs) << " MEGA k-mers per sec/proc in " << (time_second_data_pass) << " seconds" << endl;
    }

    serial_printf("%s: Total time computing load imbalance: %0.3f s\n", __FUNCTION__, time_computing_loadimbalance);
    
    CHECK_MPI(MPI_Barrier(MPI_COMM_WORLD));
    
    time_temp = MPI_Wtime();

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// GGGG: Build matrix A
/////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef NDEBUG
    {
      std::string title = "Local matrix info:";
      std::string msg   = "lrow_ids row_size: " + 
                       std::to_string(lrow_ids.size())
                       + " lcol_ids row_size: " +
                       std::to_string(lcol_ids.size())
                       + " lvals row_size: " + 
                       std::to_string(lvals.size());
      TraceUtils::print_msg(title, msg, parops);
    }
#endif

    /*! GGGG: don't think this is correct for dibella */
    assert(lrow_ids.size() == lcol_ids.size() && lcol_ids.size() == lvals.size());

    /*! Create distributed sparse matrix of sequence x kmers */
    FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, parops->grid);
    FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, parops->grid);
    FullyDistVec<uint64_t, MatrixEntry> dvals(lvals, parops->grid);

    uint64_t n_rows = seq_count;
    /*! Columns of the matrix are direct maps to kmers identified
     * by their |alphabet| base number. E.g. for proteins this is
     * base 20, so the direct map has to be 20^k in size. */

    /*! GGGG: this is the reliable k-mer space */
    auto n_cols = static_cast<uint64_t>(pow(alph.size, k));
    tp->times["start_kmerop:gen_A:spMatA()"] = std::chrono::system_clock::now();
    PSpMat<MatrixEntry>::MPI_DCCols A(n_rows, n_cols, drows, dcols, dvals, false);
    tp->times["end_kmerop:gen_A:spMatA()"] = std::chrono::system_clock::now();

    return A;
  }

  PSpMat<MatrixEntry>::MPI_DCCols KmerOps::generate_S(
      ushort k, uint64_t subk_count,
      Alphabet &alph, std::shared_ptr<ParallelOps> &parops,
      std::shared_ptr<TimePod> &tp,
      std::unordered_set<Kmer, Kmer>& local_kmers){

    dibella::Blosum62 bsm62;
    dibella::NearestKmers2 nk2(alph, bsm62);

    uvec_64 lrow_ids, lcol_ids;
    std::vector<MatrixEntry> lvals;

    tp->times["start_kmerop:gen_S:find_sub_kmers()"] = std::chrono::system_clock::now();
    for (const auto& kmer : local_kmers){
      // subk+count + 1 because we want the digonal too.
      for (size_t i = 0; i < subk_count+1; ++i){
        lrow_ids.push_back(kmer.code());
      }
      // Push diagonal column entry;
      lcol_ids.push_back(kmer.code());
      lvals.emplace_back(0, 0);

      // Now find nearest kmers.
      std::vector<Kmer> nbrs = nk2.find_sub_kmers(kmer, subk_count);
      for (auto& nbr : nbrs){
        lcol_ids.push_back(nbr.code());
        lvals.emplace_back(nbr.dist_to_root(), 0);
      }
    }
    tp->times["end_kmerop:gen_S:find_sub_kmers()"] = std::chrono::system_clock::now();
    assert(lrow_ids.size() == lcol_ids.size()
           && lcol_ids.size() == lvals.size());

    /*! Create distributed sparse matrix of sequence x kmers */
    FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, parops->grid);
    FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, parops->grid);
    FullyDistVec<uint64_t, MatrixEntry> dvals(lvals, parops->grid);

    auto n_cols = static_cast<uint64_t>(pow(alph.size, k));
    auto n_rows = n_cols;
    tp->times["start_kmerop:gen_S:spMatS()"] = std::chrono::system_clock::now();
    PSpMat<MatrixEntry>::MPI_DCCols S(n_rows, n_cols, drows, dcols, dvals, true);
    tp->times["end_kmerop:gen_S:spMatS()"] = std::chrono::system_clock::now();

    return S;
  }
}