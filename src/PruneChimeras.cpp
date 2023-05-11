#include "PruneChimeras.hpp"
#include <cassert>
#include <iostream>
#include <sstream>
#include <fstream>
#include <mpi.h>

template <class T, class T1, class T2>
T myabsdiff(T1 a, T2 b)
{
    return a < b? static_cast<T>(b - a) : static_cast<T>(b - a);
}

PileupVector::PileupVector(int read_length) : pileup(read_length, 0) {}

PileupVector::PileupVector(const std::vector<int>& v, int offset, int size) : pileup(size, 0)
{
    for (int i = 0; i < size; ++i)
        pileup[i] = v[offset + i];
}

int PileupVector::Length() const { return pileup.size(); }

void PileupVector::AddInterval(int begin, int end)
{
    assert((begin >= 0 && end <= Length()));
    for (int i = begin; i < end; ++i) pileup[i]++;
}

int PackPileupVectors(const std::vector<PileupVector>& pvs, std::vector<int>& packed, std::vector<int>& lens)
{
    int pv_size = pvs.size();
    int packed_size = 0;
    packed.clear();
    lens.reserve(pv_size);

    for (int i = 0; i < pv_size; ++i)
    {
        int len = pvs[i].Length();
        lens.push_back(len);
        packed_size += len;
    }

    packed.reserve(packed_size);

    for (int i = 0; i < pv_size; ++i)
        for (int j = 0; j < lens[i]; ++j)
            packed.push_back(pvs[i].pileup[j]);

    return packed_size;
}

void UnpackPileupVectors(const std::vector<int>& packed, const std::vector<int>& lens, std::vector<PileupVector>& pvs)
{
    int size = lens.size();
    pvs.clear();
    pvs.reserve(size);
    int offset = 0;

    for (int i = 0; i < size; ++i)
    {
        pvs.emplace_back(packed, offset, lens[i]);
        offset += lens[i];
    }
}

// void add_gaps(int begin, int end, int length, std::vector<std::tuple<int, int>>& middle, std::vector<std::tuple<int, int>>& extremity)
// {
    // if (begin != end)
    // {
        // if (!begin || end == length) extremity.push_back({begin, end});
        // else middle.push_back({begin, end});
    // }
// }

// FullyDistVec<int64_t, int64_t> GetChimeras(const std::shared_ptr<CommGrid>& commgrid, std::shared_ptr<DistributedFastaData> dfd, const std::vector<PileupVector>& pileups, int coverage_min)
// {
    // std::vector<int64_t> chimeras;
    // chimeras.reserve((pileups.size() / 10) + 1);

    // if (commgrid->GetRankInProcCol() == 0)
    // {
        // uint64_t col_seq_start_idx = dfd->col_seq_start_idx;

        // int myrank;
        // MPI_Comm_rank(commgrid->GetWorld(), &myrank);

        // for (int i = 0; i < pileups.size(); ++i)
        // {
            // const PileupVector& pv = pileups[i];

            // std::vector<std::tuple<int, int>> middle, extremity;

            // bool in_gap = true;
            // int begin = 0, end = 0;

            // for (int j = 0; j < pv.Length(); ++j)
            // {
                // if (pv.pileup[j] <= coverage_min && !in_gap)
                // {
                    // end = 0, begin = j;
                    // in_gap = true;
                // }

                // if (pv.pileup[j] > coverage_min && in_gap)
                // {
                    // end = j;
                    // in_gap = false;
                    // add_gaps(begin, end, pv.Length(), middle, extremity);
                // }
            // }

            // if (in_gap)
            // {
                // end = pv.Length();
                // add_gaps(begin, end, pv.Length(), middle, extremity);
            // }

            // if (middle.size() > 0)
            // {
                // chimeras.push_back(i+col_seq_start_idx);
                // std::cout << "Chimeric\t" << (i+col_seq_start_idx+1) << "\t";

                // for (int idx = 0; idx < middle.size(); ++idx)
                // {
                    // int g1 = std::get<0>(middle[idx]);
                    // int g2 = std::get<1>(middle[idx]);
                    // std::cout << myabsdiff<int, int, int>(g1, g2) << "," << g1 << "," << g2 << ";";
                // }
                // std::cout << std::endl;
            // }
        // }
    // }

    // return FullyDistVec<int64_t, int64_t>(chimeras, commgrid);
// }

std::vector<PileupVector> GetReadPileup(DistributedFastaData& dfd, CT<Overlap>::PSpParMat& Rmat)
{
    auto index = dfd.getindex();
    auto commgrid = index.getcommgrid();
    int myrank = commgrid->GetRank();
    int nprocs = commgrid->GetSize();

    auto spSeq = Rmat.seq();
    auto colbuf = dfd.getcolbuf();

    /* Want to calculate pileup vectors for each read. Do this by first
     * computing the pileup vectors for reads in the column range of the local
     * processor, by going through each local column and computing the pileup
     * as a function of the overlapping intervals of the column read for each
     * nonzero in that locally stored column */
    std::vector<PileupVector> local_pileups;
    size_t local_ncols = spSeq.getncol();
    assert(local_ncols == dfd.getnumcolreads());
    local_pileups.reserve(local_ncols);

    for (size_t i = 0; i < local_ncols; ++i)
    {
        int read_length = (*colbuf)[i].size();
        local_pileups.emplace_back(read_length);
    }

    size_t col_offset = dfd.getcolstartid();

    /* iterate over every local column */
    for (auto colit = spSeq.begcol(); colit != spSeq.endcol(); ++colit)
    {
        int colid = colit.colid();

        /* iterate over every row that has an overlap with current column */
        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
        {
            const Overlap& o = nzit.value();
            uint64_t rowid = nzit.rowid();
            local_pileups[colid].AddInterval(std::get<1>(o.beg), std::get<1>(o.end));
        }
    }

    std::vector<int> lens;
    std::vector<int> packed;
    MPI_Count_type size = PackPileupVectors(local_pileups, packed, lens);

    MPI_ALLREDUCE(MPI_IN_PLACE, packed.data(), size, MPI_INT, MPI_SUM, commgrid->GetColWorld());

    UnpackPileupVectors(packed, lens, local_pileups);

    return local_pileups;
}
