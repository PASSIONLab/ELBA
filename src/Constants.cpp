// Created by Saliya Ekanayake on 12/17/18.

#include "../include/Constants.hpp"
const char *CMD_OPTION_INPUT = "i";
const char *CMD_OPTION_DESCRIPTION_INPUT = "The input FASTA file.";

const char *CMD_OPTION_INPUT_SEQ_COUNT = "c";
const char *CMD_OPTION_DESCRIPTION_INPUT_SEQ_COUNT = "The sequence count in the input sequence file.";

const char *CMD_OPTION_INPUT_OVERLAP = "O";
const char *CMD_OPTION_DESCRIPTION_INPUT_OVERLAP = "Bytes to overlap when reading the input file in parallel.";

const char *CMD_OPTION_SEED_COUNT = "sc";
const char *CMD_OPTION_DESCRIPTION_SEED_COUNT = "The seed count.";

const char *CMD_OPTION_MATCH = "ma";
const char *CMD_OPTION_DESCRIPTION_MATCH = "Base match score (positive value).";

const char *CMD_OPTION_MISMATCH = "mi";
const char *CMD_OPTION_DESCRIPTION_MISMATCH = "Base mismatch penalty (negative value).";

const char *CMD_OPTION_GAP_OPEN = "g";
const char *CMD_OPTION_DESCRIPTION_GAP_OPEN = "Gap open penalty (negative value).";

const char *CMD_OPTION_GAP_EXT = "e";
const char *CMD_OPTION_DESCRIPTION_GAP_EXT = "Gap extension penalty (negative value).";

const char *CMD_OPTION_KMER_LENGTH = "k";
const char *CMD_OPTION_DESCRIPTION_KMER_LENGTH = "The k-mer length.";

const char *CMD_OPTION_KMER_STRIDE = "s";
const char *CMD_OPTION_DESCRIPTION_KMER_STRID = "The k-mer stride.";

const char *CMD_OPTION_OVERLAP_FILE = "of";
const char *CMD_OPTION_DESCRIPTION_OVERLAP_FILE = "The overlap file.";

const char *CMD_OPTION_ALIGN_FILE = "af";
const char *CMD_OPTION_DESCRIPTION_ALIGN_FILE = "The alignment file.";

const char *CMD_OPTION_NO_ALIGN = "na";
const char *CMD_OPTION_DESCRIPTION_NO_ALIGN = "Do not perform alignment.";

const char *CMD_OPTION_GPU_ALIGN = "ga";
const char *CMD_OPTION_DESCRIPTION_GPU_ALIGN = "GPU-based x-drop alignment.";

const char *CMD_OPTION_CPU_ALIGN = "ca";
const char *CMD_OPTION_DESCRIPTION_CPU_ALIGN = "CPU-based x-drop alignment.";

const char *CMD_OPTION_IDX_MAP = "idxmap";
const char *CMD_OPTION_DESCRIPTION_IDX_MAP = "The file path to write the global sequence indices to original global sequence indices mapping.";

const char *CMD_OPTION_ALPH = "alph";
const char *CMD_OPTION_DESCRIPTION_ALPH = "The alphabet to use - valid choices are [dna|protein].";

const char *CMD_OPTION_SUBS = "subs";
const char *CMD_OPTION_DESCRIPTION_SUBS = "# substitute kmers.";

const char *CMD_OPTION_JOB_NAME_PREFIX = "jp";
const char *CMD_OPTION_DESCRIPTION_JOB_NAME_PREFIX = "Job name prefix.";

const char *CMD_OPTION_LOG_FREQ = "lf";
const char *CMD_OPTION_DESCRIPTION_LOG_FREQ = "Log frequency.";

const char *CMD_OPTION_AF_FREQ = "afreq";
const char *CMD_OPTION_DESCRIPTION_AF_FREQ = "The alignment write frequency in number of lines.";