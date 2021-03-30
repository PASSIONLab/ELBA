#include <cmath>
#include "../include/Constants.hpp"
#include "../include/ParallelOps.hpp"
#include "../include/ParallelFastaReader.hpp"
#include "../include/Alphabet.hpp"
#include "../include/DistributedPairwiseRunner.hpp"
#include "../include/cxxopts.hpp"
#include "../include/TransitiveReductionSR.hpp"
#include "../include/Utils.hpp"

#include <map>
#include <fstream>

/*! Namespace declarations */
using namespace combblas;

void GetAssembly(std::vector<std::string>& myContigSet, TraceUtils tu)
{
    /* @GGGG-TODO: Contig correction and scaffolding() 
     * (1) Partial order alignment of excluded sequences (that passed aligment threshold) to contig to correct 
     * them (keep track haplotypes that are not errors based on the coverage)
     * (2) |contig|x|contig| jaccard similarity (or some sort of similarity) to scaffold them 
     * */
}