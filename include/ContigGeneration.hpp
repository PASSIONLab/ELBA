#ifndef CONTIG_GENERATION_HPP_
#define CONTIG_GENERATION_HPP_

#include <vector>
#include <string>
#include "Overlap.hpp"
#include "DistributedFastaData.hpp"
#include "common.h"

std::vector<std::string> GenerateContigs(CT<Overlap>::PSpParMat& S, const DnaBuffer& mydna, DistributedFastaData& dfd);

#endif
