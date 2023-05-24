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

std::tuple<int, int> PileupVector::GetTrimmedInterval(int threshold)
{
    int len = Length();
    int beststart, bestend;
    double bestavg = 0;
    int curbases = 0;
    int start = -1, end = -1, maxlen = 2500;

    for (int i = 0; i < len; ++i)
    {
        if (pileup[i] >= threshold)
        {
            if (start == -1)
            {
                curbases = 0;
                start = i;
            }

            end = i;
            curbases += pileup[i];
            int span = end - start + 1;
            double curavg = static_cast<double>(curbases) / static_cast<double>(span);

            if (span > maxlen && curavg > bestavg)
            {
                beststart = start;
                bestend = end;
                maxlen = span;
                bestavg = curavg;
            }
        }
        else
        {
            start = -1;
            end = -1;
        }
    }

    return {start, end};
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
        int64_t colid = colit.colid();

        /* iterate over every row that has an overlap with current column */
        for (auto nzit = spSeq.begnz(colit); nzit != spSeq.endnz(colit); ++nzit)
        {
            const Overlap& o = nzit.value();
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


