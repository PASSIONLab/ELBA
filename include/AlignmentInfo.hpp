// Created by Saliya Ekanayake on 2019-07-10.

#ifndef LBL_PISA_ALIGNMENTINFO_HPP
#define LBL_PISA_ALIGNMENTINFO_HPP

struct AlignmentInfo{
  seqan::AlignmentStats stats;
  ushort seq_h_length;
  ushort seq_v_length;
  ushort seq_h_seed_length;
  ushort seq_v_seed_length;
  uint64_t seq_h_g_idx;
  uint64_t seq_v_g_idx;
};

#endif //LBL_PISA_ALIGNMENTINFO_HPP
