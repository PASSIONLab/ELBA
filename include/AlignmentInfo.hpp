// Created by Saliya Ekanayake on 2019-07-10 and modified by Giulia Guidi on 08/31/2020.

#ifndef ELBA_ALIGNMENTINFO_HPP
#define ELBA_ALIGNMENTINFO_HPP

#include <seqan/align.h>
#include "Types.hpp"

struct AlignmentInfo
{

  // For others, not x-drop
  seqan::AlignmentStats stats;

	int xscore;
	TSeed seed;
  
  bool rc;

  ushort seq_h_length;
  ushort seq_v_length;
  
  ushort seq_h_seed_length;
  ushort seq_v_seed_length;

  uint64_t seq_h_g_idx;
  uint64_t seq_v_g_idx;
};

#endif //ELBA_ALIGNMENTINFO_HPP
