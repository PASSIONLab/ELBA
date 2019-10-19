// Created by Saliya Ekanayake on 10/15/19.

#include "../../include/kmer/KmerOps.hpp"
#include "../../include/NearestKmers2.hpp"

namespace pisa {
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

    tp->times["start_kmerop:gen_A:loop_add_kmers()"] = std::chrono::system_clock::now();
    for (uint64_t lseq_idx = 0; lseq_idx < lfd->local_count(); ++lseq_idx) {
      buff = lfd->get_sequence(lseq_idx, len, start_offset,
                               end_offset_inclusive);
      auto num_kmers = add_kmers(buff, len, start_offset, end_offset_inclusive,
                                 k, s, alph, lcol_ids, lvals, parops, local_kmers);
      lrow_ids.insert(lrow_ids.end(), num_kmers, lseq_idx + offset);
    }
    tp->times["end_kmerop:gen_A:loop_add_kmers()"] = std::chrono::system_clock::now();


#ifndef NDEBUG
    {
      std::string title = "Local matrix info:";
      std::string msg = "lrow_ids row_size: " + std::to_string(lrow_ids.size())
                        + " lcol_ids row_size: " +
                        std::to_string(lcol_ids.size())
                        + " lvals row_size: " + std::to_string(lvals.size());
      TraceUtils::print_msg(title, msg, parops);
    }
#endif

    assert(lrow_ids.size() == lcol_ids.size()
    && lcol_ids.size() == lvals.size());

    /*! Create distributed sparse matrix of sequence x kmers */
    FullyDistVec<uint64_t, uint64_t> drows(lrow_ids, parops->grid);
    FullyDistVec<uint64_t, uint64_t> dcols(lcol_ids, parops->grid);
    FullyDistVec<uint64_t, MatrixEntry> dvals(lvals, parops->grid);

    uint64_t n_rows = seq_count;
    /*! Columns of the matrix are direct maps to kmers identified
     * by their |alphabet| base number. E.g. for proteins this is
     * base 20, so the direct map has to be 20^k in size. */

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

    pisa::Blosum62 bsm62;
    pisa::NearestKmers2 nk2(alph, bsm62);

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