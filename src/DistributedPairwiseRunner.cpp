//
// Created by Saliya Ekanayake on 2019-02-19.
// Modified by Aydin Buluc on 2019-12-29 
//


#include "../include/DistributedPairwiseRunner.hpp"
#include <atomic>         // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT


DistributedPairwiseRunner::DistributedPairwiseRunner(
    const std::shared_ptr<DistributedFastaData> dfd,
    PSpMat<distal::CommonKmers>::DCCols * localmat,
    int afreq,
    uint64_t rowoffset, uint64_t coloffset,
    const std::shared_ptr<ParallelOps> &parops)
    : dfd(dfd), spSeq(localmat), row_offset(rowoffset), col_offset(coloffset), afreq(afreq), parops(parops) {

}


void DistributedPairwiseRunner::write_overlaps(const char *file) {

  uint64_t local_nnz_count = 0;
  uint64_t local_top_triangle_count = 0;
  std::stringstream ss;

  if (parops->world_proc_rank == 0) {
    ss << "g_col_idx,g_row_idx,common_kmer_count" << std::endl;
  }
  ushort l_max_common_kmers = 0;
  for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit) {
    // iterate over columns
    auto l_col_idx = colit.colid(); // local numbering
    uint64_t g_col_idx = l_col_idx + col_offset;

    for (auto nzit = spSeq->begnz(colit); nzit < spSeq->endnz(colit); ++nzit) {
      auto l_row_idx = nzit.rowid();
      uint64_t g_row_idx = l_row_idx + row_offset;

      ++local_nnz_count;

      /*!
       * Note. the cells means the process grid cells.
       * We only want to compute the top triangles of any grid cell.
       * Further, we want cell diagonals only for cells that are on the
       * top half of the grid excluding the grid's main diagonal cells
       */
      if (l_col_idx < l_row_idx){
        continue;
      }

      if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx){
        continue;
      }

      distal::CommonKmers cks = nzit.value();
      if (cks.count > l_max_common_kmers){
        l_max_common_kmers = cks.count;
      }

      ++local_top_triangle_count;
      ss << g_col_idx << "," << g_row_idx << "," << cks.count << "\n";
    }
  }

  ushort g_max_common_kmers = 0;
  MPI_Reduce(&l_max_common_kmers, &g_max_common_kmers, 1,
             MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);
  if (parops->world_proc_rank == 0){
    std::printf("  Max common kmers %d\n", g_max_common_kmers);
  }
  std::string overlaps_str = ss.str();
  parops->write_file_in_parallel(file, overlaps_str);
}

void DistributedPairwiseRunner::run(PairwiseFunction *pf, const char* file, std::ofstream& lfs, int log_freq) {
  /*! There are two types of rows and columns below.
   * The sequences are arranged as an NxN matrix in
   * mat (this is not how it's stored internally).
   * This NxN is distributed over a grid of
   * sqrt(P) x sqrt (P), where P is the total number
   * of processes. Anything to do with the grid will
   * be prefixed by gr_*/


  std::ofstream af_stream;
  af_stream.open(file);

  uint64_t local_nnz_count = spSeq->getnnz();
  std::atomic<uint64_t> current_nnz_count(0);

  lfs << "Local nnz count: " << local_nnz_count << std::endl;

  int numThreads = 1;	// default case
#ifdef THREADED
#pragma omp parallel
    {
        numThreads = omp_get_num_threads();
    }
#endif

  std::vector<std::stringstream> ss(numThreads);
  if(parops->world_proc_rank == 0){
    af_stream << "g_col_idx,g_row_idx,pid,col_seq_len,row_seq_len,col_seq_align_len,row_seq_align_len, num_gap_opens, col_seq_len_coverage, row_seq_len_coverage" << std::endl;
  }

  std::atomic<uint64_t> line_count(0);

  PSpMat<pisa::CommonKmers>::Tuples mattuples(*spSeq);
 
#pragma omp parallel for
  for(uint64_t i=0; i< local_nnz_count; i++)
  {
	  auto l_row_idx = mattuples.rowindex(i);
	  auto l_col_idx = mattuples.colindex(i);
    	  uint64_t g_col_idx = l_col_idx + col_offset;
	  uint64_t g_row_idx = l_row_idx + row_offset;		  

	  seqan::Peptide *seq_h = dfd->col_seq(l_col_idx);  
	  seqan::Peptide *seq_v = dfd->row_seq(l_row_idx);

	  current_nnz_count++;
	  if (current_nnz_count % log_freq == 0){
		  #pragma omp critical
		  {
        		auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
        		lfs << "  (" << current_nnz_count << "/" << local_nnz_count << ") -- "
        		<< std::setprecision(2) << (1.0*current_nnz_count / local_nnz_count)
        		<< "% done. " << std::ctime(&t);
        		lfs.flush();
		  }
	  }	

	  /*!
	   * Note. the cells means the process grid cells.
	   * We only want to compute the top triangles of any grid cell.
	   * Further, we want cell diagonals only for cells that are on the
	   * top half of the grid excluding the grid's main diagonal cells
	   * */
	  if (l_col_idx < l_row_idx){
		  continue;
	  }
	  if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx){
		  continue;
	  }

	  pisa::CommonKmers cks = mattuples.numvalue(i);


	  int myThread = 0;
#ifdef THREADED
	  myThread = omp_get_thread_num();
#endif

	  pf->apply(l_col_idx, g_col_idx, l_row_idx, g_row_idx, seq_h, seq_v, cks, ss[myThread]);
	  line_count++;

	  if (line_count%afreq == 0){
		  #pragma omp critical
		  {
		  	af_stream << ss[myThread].str();
		  	af_stream.flush();
		  	ss[myThread].str(std::string());
		  }
	  }
  }

  lfs << "  (" << current_nnz_count << "/" << local_nnz_count << ") -- "
      << "100% done." << std::endl;

  pf->print_avg_times(parops);
//  if(parops->world_proc_rank == 0) {
//  }
//  ushort g_max_common_kmers = 0;
//  MPI_Reduce(&l_max_common_kmers, &g_max_common_kmers, 1,
//             MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);
//  if (parops->world_proc_rank == 0){
//    std::printf("  Max common kmers %d\n", g_max_common_kmers);
//  }
//  std::string align_str = ss.str();
//  parops->write_file_in_parallel(file, align_str);

  for(int i=0; i< numThreads; ++i)
  {
  	af_stream << ss[i].str();
  }
  af_stream.flush();
  af_stream.close();
}
