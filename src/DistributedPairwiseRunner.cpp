//
// Created by Saliya Ekanayake on 2019-02-19
// Modified by Aydin Buluc on 2019-12-29 
// Modified by Gilia Guidi on 2021-03-09 
//

#include "../include/DistributedPairwiseRunner.hpp"
#include <atomic> // std::atomic, std::atomic_flag, ATOMIC_FLAG_INIT

DistributedPairwiseRunner::DistributedPairwiseRunner(
    const std::shared_ptr<DistributedFastaData> dfd,
    PSpMat<elba::CommonKmers>::DCCols * localmat,
	PSpMat<elba::CommonKmers>::MPI_DCCols * glmat,
    int afreq,
    uint64_t rowoffset, uint64_t coloffset,
    const std::shared_ptr<ParallelOps> &parops)
    : dfd(dfd), gmat(glmat), spSeq(localmat), row_offset(rowoffset),
			col_offset(coloffset), afreq(afreq), parops(parops) {
}

void DistributedPairwiseRunner::write_overlaps(const char *file)
{
	uint64_t local_nnz_count = 0;
	uint64_t local_top_triangle_count = 0;
	std::stringstream ss;

	std::ofstream afs;
	afs.open(file);

	if (parops->world_proc_rank == 0)
	{
	  ss << "g_col_idx,g_row_idx,common_kmer_count" << std::endl;
	}

	ushort l_max_common_kmers = 0;
	for (auto colit = spSeq->begcol(); colit != spSeq->endcol(); ++colit)
	{
	  // iterate over columns
	  auto l_col_idx = colit.colid(); // local numbering
	  uint64_t g_col_idx = l_col_idx + col_offset;

	  for (auto nzit = spSeq->begnz(colit); nzit < spSeq->endnz(colit); ++nzit)
		{
	    auto l_row_idx = nzit.rowid();
	    uint64_t g_row_idx = l_row_idx + row_offset;

	    ++local_nnz_count;

	    /*!
	     * Note. the cells means the process grid cells.
	     * We only want to compute the top triangles of any grid cell.
	     * Further, we want cell diagonals only for cells that are on the
	     * top half of the grid excluding the grid's main diagonal cells
	     */
	    if (l_col_idx < l_row_idx)
		  {
	      continue;
	    }

	    if (l_col_idx == l_row_idx && g_col_idx <= g_row_idx)
		  {
	      continue;
	    }

	    elba::CommonKmers cks = nzit.value();
	    if (cks.count > l_max_common_kmers)
		{
	      l_max_common_kmers = cks.count;
	    }

	    ++local_top_triangle_count;
		  afs << g_row_idx << " " << g_col_idx << "\n";
	  }
	}

	afs.close();

	ushort g_max_common_kmers = 0;
	MPI_Reduce(&l_max_common_kmers, &g_max_common_kmers, 1,
	           MPI_UINT16_T, MPI_MAX, 0, MPI_COMM_WORLD);

	if (parops->world_proc_rank == 0)
	{
	  std::printf("  max common kmers %d\n", g_max_common_kmers);
	}
}

void DistributedPairwiseRunner::run(PairwiseFunction *pf, const char* file, std::ofstream& lfs, int log_freq, ushort k) {
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
  if(parops->world_proc_rank == 0)
    af_stream << "g_col_idx,g_row_idx,pid,col_seq_len,row_seq_len,"
		"col_seq_align_len,row_seq_align_len, num_gap_opens, "
		"col_seq_len_coverage, row_seq_len_coverage,common_count" << std::endl;

  std::atomic<uint64_t> line_count(0);
  uint64_t nalignments = 0;
  PSpMat<elba::CommonKmers>::Tuples mattuples(*spSeq);
  
  #pragma omp parallel for reduction(+:nalignments)
  for(uint64_t i=0; i< local_nnz_count; i++)
  {
	auto l_row_idx = mattuples.rowindex(i);
	auto l_col_idx = mattuples.colindex(i);
	uint64_t g_col_idx = l_col_idx + col_offset;
	uint64_t g_row_idx = l_row_idx + row_offset;		  

	seqan::Dna5String *seq_h = dfd->col_seq(l_col_idx);  
	seqan::Dna5String *seq_v = dfd->row_seq(l_row_idx);

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

	elba::CommonKmers cks = mattuples.numvalue(i);

	int myThread = 0;
#ifdef THREADED
	myThread = omp_get_thread_num();
#endif

		++nalignments;
		pf->apply(l_col_idx, g_col_idx, l_row_idx, g_row_idx, seq_h, seq_v, k, cks, ss[myThread]);
		line_count++;

		if (line_count % afreq == 0)
		{
			#pragma omp critical
			{
			  af_stream << ss[myThread].str();
			  af_stream.flush();
			  ss[myThread].str(std::string());
			}
		}
	}

	pf->nalignments = nalignments;

	lfs << "  (" << current_nnz_count << "/" << local_nnz_count << ") -- "
	    << "100% done." << std::endl;
	lfs << "#alignments run " << nalignments << std::endl;

	pf->print_avg_times(parops, lfs);
	
	for(int i=0; i< numThreads; ++i)
	{
		af_stream << ss[i].str();
	}
	af_stream.flush();
	af_stream.close();
}


void
DistributedPairwiseRunner::run_batch
(
    PairwiseFunction	*pf,
	std::ofstream&		 lfs,
	int					 log_freq,
	int					 ckthr,
	bool				 aln_score_thr,
	TraceUtils 			 tu,
	const bool 			 noAlign,
	ushort 				 k,
	uint64_t 			 nreads,
	bool				 score_only
)
{
	uint64_t	local_nnz_count = spSeq->getnnz();
	
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int			batch_size		= 1e5;
	int			batch_cnt		= (local_nnz_count / batch_size) + 1;
	int			batch_idx		= 0;
	uint64_t	nalignments		= 0;

	// PSpMat<elba::CommonKmers>::Tuples mattuples(*spSeq);
	// @TODO threaded
	PSpMat<elba::CommonKmers>::ref_tuples *mattuples =
		new PSpMat<elba::CommonKmers>::ref_tuples[local_nnz_count];

	uint64_t z = 0;
	auto dcsc = spSeq->GetDCSC();

	for (uint64_t i = 0; i < dcsc->nzc; ++i)
	{
		for (uint64_t j = dcsc->cp[i]; j < dcsc->cp[i+1]; ++j)
		{
			std::get<0>(mattuples[z]) = dcsc->ir[j];
			std::get<1>(mattuples[z]) = dcsc->jc[i];
			std::get<2>(mattuples[z]) = &(dcsc->numx[j]);
			++z;
		}
	}

	assert (z == local_nnz_count);
		
	lfs << "Local nnz count: " << local_nnz_count << std::endl;

	int numThreads = 1;
	#ifdef THREADED
	#pragma omp parallel
    {
      	numThreads = omp_get_num_threads();
    }
	#endif

	uint64_t *algn_cnts   = new uint64_t[numThreads + 1];
	uint64_t nelims_ckthr = 0; // nelims_alnthr = 0, nelims_both = 0;

	// GGGG: local vector containing global indexes of sequences whose nonzeros should be removed because contained vertex
	std::vector<std::vector<int64_t>> ContainedSeqPerBatch(batch_cnt);
	
	while (batch_idx < batch_cnt) 
	{
		uint64_t beg = batch_idx * batch_size;
		uint64_t end = ((batch_idx + 1) * batch_size > local_nnz_count) ? local_nnz_count : ((batch_idx + 1) * batch_size);

	#ifdef VERBOSE
		tu.print_str("Batch idx " + std::to_string(batch_idx) + "/" +
					 std::to_string(batch_cnt) + " [" +
					 std::to_string(beg) + ", " +
					 std::to_string(end) + ")\n");
	#endif

		memset(algn_cnts, 0, sizeof(*algn_cnts) * (numThreads + 1));

		uint64_t nelims_ckthr_cur = 0;
		
		// Count number of alignments in this batch
		#pragma omp parallel reduction(+:nelims_ckthr_cur)
		{
			int tid = 0;
			#ifdef THREADED
			tid = omp_get_thread_num();
			#endif

			uint64_t algn_cnt = 0;

			#pragma omp for schedule(static, 1000)
			for (uint64_t i = beg; i < end; ++i)
			{
				auto				 l_row_idx = std::get<0>(mattuples[i]);
				auto				 l_col_idx = std::get<1>(mattuples[i]);
				uint64_t			 g_col_idx = l_col_idx + col_offset;
				uint64_t			 g_row_idx = l_row_idx + row_offset;

				assert(l_row_idx >= 0 && l_col_idx >= 0 && g_col_idx >= 0 && g_row_idx >= 0);

				elba::CommonKmers *cks = std::get<2>(mattuples[i]);

				if ((cks->count >= ckthr) 	 	&& 
					(l_col_idx >= l_row_idx) 	&&
					(l_col_idx != l_row_idx  || g_col_idx > g_row_idx))
				{
					++algn_cnt;
				}

				if ((l_col_idx >= l_row_idx) && (l_col_idx != l_row_idx || g_col_idx > g_row_idx))
				{
					if (cks->count < ckthr) ++nelims_ckthr_cur;
				}
			}

			algn_cnts[tid + 1] = algn_cnt;
		}

		nelims_ckthr += nelims_ckthr_cur;	

		for (int i = 1; i < numThreads + 1; ++i) algn_cnts[i] += algn_cnts[i - 1];

		nalignments += algn_cnts[numThreads];

		if (algn_cnts[numThreads] == 0)
		{	
			++batch_idx;
			continue;
		}
		
		// allocate StringSet
		seqan::StringSet<seqan::Dna5String> seqsh;
		seqan::StringSet<seqan::Dna5String> seqsv;
		resize(seqsh, algn_cnts[numThreads], seqan::Exact{});
		resize(seqsv, algn_cnts[numThreads], seqan::Exact{});

		uint64_t *lids = new uint64_t[algn_cnts[numThreads]];
		
		// fill StringSet
		#pragma omp parallel
		{
			int tid = 0;
			#ifdef THREADED
			tid = omp_get_thread_num();
			#endif

			uint64_t algn_idx = algn_cnts[tid];

			#pragma omp for schedule(static, 1000)
			for (uint64_t i = beg; i < end; ++i)
			{
				auto		l_row_idx = std::get<0>(mattuples[i]);
				auto		l_col_idx = std::get<1>(mattuples[i]);
				uint64_t	g_col_idx = l_col_idx + col_offset;
				uint64_t	g_row_idx = l_row_idx + row_offset; 

				assert(l_row_idx >= 0 && l_col_idx >= 0 && g_col_idx >= 0 && g_row_idx >= 0);

				elba::CommonKmers *cks = std::get<2>(mattuples[i]);

				if ((cks->count >= ckthr) && (l_col_idx >= l_row_idx) && (l_col_idx != l_row_idx  || g_col_idx > g_row_idx))
				{

					seqsh[algn_idx] = seqan::Dna5String(*(dfd->col_seq(l_col_idx)));
					seqsv[algn_idx] = seqan::Dna5String(*(dfd->row_seq(l_row_idx)));

					lids[algn_idx] = i;
					++algn_idx;
				}
			}
		}

		// function call to the aligner
		lfs << "calling aligner for batch idx " << batch_idx
			<< " cur #algnments " 				<< algn_cnts[numThreads]
			<< " overall " 						<< nalignments
			<< std::endl;

		// GGGG: fill ContainedSeqPerBatch
		pf->apply_batch(seqsh, seqsv, lids, col_offset, row_offset, mattuples, lfs, noAlign, k, nreads, ContainedSeqPerBatch[batch_idx]);
		
		delete [] lids;
		++batch_idx;
	}

	// lfs << "Before concatenation ContainedSeqPerBatch" << std::endl;
	
	int readcount = 0;
	for(int b = 0; b < batch_cnt; ++b)
	{
		readcount += ContainedSeqPerBatch[b].size();
	}

	unsigned int readssofar = 0; 
	std::vector<int64_t> ContainedSeqPerProc(readcount);

	// Concatenate per-batch result
	for(int b = 0; b < batch_cnt; ++b)
	{
		copy(ContainedSeqPerBatch[b].begin(), ContainedSeqPerBatch[b].end(), ContainedSeqPerProc.begin() + readssofar);
		readssofar += ContainedSeqPerBatch[b].size();
	}

	lfs << "Post concatenation ContainedSeqPerBatch, readcount " << readcount << std::endl;

	std::sort(ContainedSeqPerProc.begin(), ContainedSeqPerProc.end());
	ContainedSeqPerProc.erase(unique(ContainedSeqPerProc.begin(), ContainedSeqPerProc.end()), ContainedSeqPerProc.end());

	lfs << "Removed duplicates and sorted, contained readcount " << ContainedSeqPerProc.size() << " out of " << nreads << " sequences" << std::endl;

	//////////////////////////////////////////////////////////////////////////////////////
	// CONTAINED SEQRUENCES COMMUNICATION                                               // 
	//////////////////////////////////////////////////////////////////////////////////////

	// Don't use boolean (it's bitmap not array of boolean, this messes up in communication)
	FullyDistVec<int64_t, int64_t> ContainedSeqGlobal(gmat->getcommgrid(), nreads, 0); // I need to use same type as index for Prune()

	FullyDistVec<int64_t, int64_t> ri(gmat->getcommgrid(), nreads, 0); // I need to use same type as index for Prune()
	FullyDistVec<int64_t, int64_t> ci(gmat->getcommgrid(), nreads, 0); // I need to use same type as index for Prune()

	int nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// A vector of vector for communication
	std::vector<std::vector<int64_t>> buffer(nprocs); // outer dim is the number of processes
	
	int * sendcnt = new int[nprocs](); // zero initialize
	int * recvcnt = new int[nprocs];

	// Each proc has a bunch of int64_t vectors with indexes (possible duplicates intra- and inter-proc)
	for(int i = 0; i < ContainedSeqPerProc.size(); i++)
	{
		int64_t lid; 

		//! Given global index gind,
		//! Return the owner processor id, and
		//! Assign the local index to lind
		int owner = ContainedSeqGlobal.Owner(ContainedSeqPerProc[i], lid);  // find the owner (global_index, local_index)

		int64_t gid = ContainedSeqPerProc[i];

		buffer[owner].push_back(gid);
		++sendcnt[owner];
	}

	// GGGG: Alltoall to share the count and Alltoallv for vector 
	MPI_Alltoall(sendcnt, 1, MPI_INT, recvcnt, 1, MPI_INT, MPI_COMM_WORLD);

	std::vector<int64_t> sendbuffer;
	std::vector<int64_t> recvbuffer;

	// Concatenate buffer for Alltoallv
	for(int i = 0; i < nprocs; ++i)
	{
		sendbuffer.insert(sendbuffer.end(), buffer[i].begin(), buffer[i].end());
	}

	// Compute displacements for Alltoallv
	int * sdispls = new int[nprocs];
	int * rdispls = new int[nprocs];

	sdispls[0] = 0;
	rdispls[0] = 0;

	for(int i = 0; i < nprocs-1; ++i)
	{
		sdispls[i+1] = sdispls[i] + sendcnt[i];
		rdispls[i+1] = rdispls[i] + recvcnt[i];
	}

	int64_t totrecv = std::accumulate(recvcnt, recvcnt + nprocs, static_cast<int64_t>(0));
	recvbuffer.resize(totrecv);
	
	// GGGG: Alltoallv for vector
	// Data is already in the right order (recepeint will receive from nrpocs processes)
	MPI_Alltoallv(sendbuffer.data(), sendcnt, sdispls, MPI_LONG, recvbuffer.data(), recvcnt, rdispls, MPI_LONG, MPI_COMM_WORLD);
	DeleteAll(sendcnt, recvcnt, sdispls, rdispls);
	
	// Go over recvbuffer and set elements to 1 in ContainedSeqPerGlobal
	for(int k = 0; k < totrecv; k++)
	{
		int idx = recvbuffer[k];
		ContainedSeqGlobal.SetElement(idx, 1); // SetElement only work locally (owner)
	}

	// Use FindInds https://github.com/PASSIONLab/CombBLAS/blob/master/include/CombBLAS/FullyDistSpVec.cpp
	// Pass the returned vec and then pass it to the Prune function (two identical vec)
	FullyDistVec<int64_t, int64_t> toerase = ContainedSeqGlobal.FindInds(bind2nd(greater<int64_t>(), 0));	// Only the non-zero indices (contained sequences)
	// toerase.DebugPrint(); // 0 is here but it's not erase from the matrix

	//////////////////////////////////////////////////////////////////////////////////////
	// TOTAL ALIGNMENTES COMMUNICATION                                                  // 
	//////////////////////////////////////////////////////////////////////////////////////

	if(noAlign) nalignments = 0;

	pf->nalignments = nalignments;
	pf->print_avg_times(parops, lfs);
 
	lfs << "#alignments run " << nalignments << std::endl;

	// Compute statistics
	uint64_t nelims_ckthr_tot = 0, nalignments_tot = 0, maxalignments = 0, minalignments = 0;

	MPI_Reduce(&nelims_ckthr, &nelims_ckthr_tot, 1, MPI_UINT64_T,
			   MPI_SUM, 0, MPI_COMM_WORLD);

	MPI_Reduce(&nalignments, &nalignments_tot, 1, MPI_UINT64_T,
			   MPI_SUM, 0, MPI_COMM_WORLD);


	uint64_t avgalignments = nalignments_tot / parops->world_procs_count;

	// min, max num alignments per proc
  	MPI_Reduce(&nalignments, &maxalignments, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_WORLD);
 	MPI_Reduce(&nalignments, &minalignments, 1, MPI_UINT64_T, MPI_MIN, 0, MPI_COMM_WORLD);

	tu.print_str("#nonzeros in C before pruning: ");
    tu.print_str(std::to_string(gmat->getnnz()) + "\n");

	//////////////////////////////////////////////////////////////////////////////////////
	// PRUNE MATRIX FROM SPURIOUS AND CONTAINED ALIGNMENT                               // 
	//////////////////////////////////////////////////////////////////////////////////////

	// Prune entries that do not meet score criteria
	auto elim_score = [] (elba::CommonKmers &ck) { return ck.passed == false; };
	gmat->Prune(elim_score); 

	// GGGG: if noAlign == true, we remove only the contained overlaps as they are not useful for transitive reduction (next prune)
	tu.print_str("#nonzeros in C after score pruning: " + std::to_string(gmat->getnnz()) + "\n");

#ifdef VERBOSE
    toerase.PrintInfo("toerase");
#endif 

	// Prune entries involving contained sequences
	gmat->PruneFull(toerase, toerase);

	tu.print_str("#nonzeros in C after contained sequences pruning: " + std::to_string(gmat->getnnz()) + "\n");
	
	delete [] algn_cnts;
	delete [] mattuples;

	return;
}
