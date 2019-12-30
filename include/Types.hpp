// Created by Saliya Ekanayake on 2/7/19.

#ifndef LBL_DAL_TYPES_HPP
#define LBL_DAL_TYPES_HPP

#include <iostream>
#include <vector>
#include <ratio>
#include <chrono>
#include <seqan/align.h>
#include <seqan/seeds.h>

typedef std::chrono::duration<double, std::milli> ms_t;
typedef std::chrono::duration<double, std::micro> micros_t;
typedef std::chrono::time_point<std::chrono::high_resolution_clock> ticks_t;
typedef std::chrono::high_resolution_clock hrc_t;

typedef seqan::Seed<seqan::Simple> TSeed;

typedef unsigned short ushort;
typedef unsigned char uchar;
typedef std::vector<uint64_t> uvec_64;
/*! TODO - apparently there's a bug when setting a different element type,
 * so let's use uint64_t as element type for now
 */
//typedef std::vector<ushort> uvec_16;
typedef std::vector<uint64_t> uvec_16;

#endif //LBL_DAL_TYPES_HPP
