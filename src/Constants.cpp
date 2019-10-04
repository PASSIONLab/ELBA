// Created by Saliya Ekanayake on 12/17/18.

#include "../include/Constants.hpp"
const char *CMD_OPTION_INPUT = "i";
const char *CMD_OPTION_DESCRIPTION_INPUT = "The input FASTA file.";

const char *CMD_OPTION_INPUT_SEQ_COUNT = "c";
const char *CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT = "The sequence count in the input sequence file.";

const char *CMD_OPTION_INPUT_OVERLAP = "O";
const char *CMD_OPTION_DESCRIPTION_INPUT_OVERLAP = "Number of bytes to overlap when reading the input file in parallel.";

const char *CMD_OPTION_SEED_COUNT = "sc";
const char *CMD_OPTION_DESCRIPTION_SEED_COUNT = "Seed count.";

const char *CMD_OPTION_GAP_OPEN = "g";
const char *CMD_OPTION_DESCRIPTION_GAP_OPEN = "Gap open penalty (negative value).";

const char *CMD_OPTION_GAP_EXT = "e";
const char *CMD_OPTION_DESCRIPTION_GAP_EXT = "Gap extension penalty (negative value).";

const char *CMD_OPTION_KMER_LENGTH = "k";
const char *CMD_OPTION_DESCRIPTION_KMER_LENGTH = "Kmer length.";

const char *CMD_OPTION_KMER_STRIDE = "s";
const char *CMD_OPTION_DESCRIPTION_KMER_STRID = "Kmer stride.";

const char *CMD_OPTION_OVERLAP_FILE = "of";
const char *CMD_OPTION_DESCRIPTION_OVERLAP_FILE = "Overlap file.";

const char *CMD_OPTION_ALIGN_FILE = "af";
const char *CMD_OPTION_DESCRIPTION_ALIGN_FILE = "Alignment file.";

const char *CMD_OPTION_NO_ALIGN = "na";
const char *CMD_OPTION_DESCRIPTION_NO_ALIGN = "Flag to indicate not to performa alignments.";

const char *CMD_OPTION_FULL_ALIGN = "fa";
const char *CMD_OPTION_DESCRIPTION_FULL_ALIGN = "Flag to indicate full alignment";

const char *CMD_OPTION_XDROP_ALIGN = "xa";
const char *CMD_OPTION_DESCRIPTION_XDROP_ALIGN = "Flag to indicate seed-and-extend using xdrop alignment";

const char *CMD_OPTION_BANDED_ALIGN = "ba";
const char *CMD_OPTION_DESCRIPTION_BANDED_ALIGN = "Flag to indicate banded alignment";

const char *CMD_OPTION_IDX_MAP = "idxmap";
const char *CMD_OPTION_DESCRIPTION_IDX_MAP = "The file path to write the global sequence indices to original global sequence indices mapping. ";

const char *CMD_OPTION_ALPH = "alph";
const char *CMD_OPTION_DESCRIPTION_ALPH = "The alphabet to use. Valid choices are [protein].";

const char *CMD_OPTION_SUBS = "subs";
const char *CMD_OPTION_DESCRIPTION_SUBS = "Add substitute kmers.";

const char *CMD_OPTION_JOB_NAME_PREFIX = "jp";
const char *CMD_OPTION_DESCRIPTION_JOB_NAME_PREFIX = "Job name prefix.";

const char *CMD_OPTION_LOG_FREQ = "lf";
const char *CMD_OPTION_DESCRIPTION_LOG_FREQ = "Log frequency.";