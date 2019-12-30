// Created by Saliya Ekanayake on 11/26/19.
#include <iostream>
#include <string>
#include <seqan/align.h>

typedef std::chrono::duration<double, std::milli> ms_t;

int main(int argc, char** argv){
  std::string str1("mvlsegewqlvlhvwakveadvaghgqdilirlfkshpetlekfdrvkhlkteaemkasedlkkhgvtvltalgailkkkghheaelkplaqshatkhkipikylefiseaiihvlhsrhpgnfgadaqgamnkalelfrkdiaakykelgyqg");
  std::string str2("mvlsegewqlvlhvwakveadvaghgqdilirlfkshpetlekfdrfkhlkteaemkasedlkkagvtvltalgailkkkghheaelkplaqshatkhkipikylefiseaiihvlhsrhpgnfgadaqgamnkalelfrkdiaakykelgyqg");

  seqan::Peptide seq1(str1);
  seqan::Peptide seq2(str2);

  std::cout << seq1 << std::endl;

  int gap_ext = -2;
  int gap_open = -11;
  seqan::Blosum62 blosum62(gap_ext, gap_open);

  seqan::Align<seqan::Peptide> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq1);

  auto start_time = std::chrono::system_clock::now();
  int score = localAlignment(align, blosum62, seqan::DynamicGaps());
  auto end_time = std::chrono::system_clock::now();
  std::cout << "FA:local_alignment" <<  (ms_t(end_time - start_time)).count() << "ms" << std::endl;

  start_time = std::chrono::system_clock::now();
  seqan::AlignmentStats stats;
  computeAlignmentStats(stats, align, blosum62);
  end_time = std::chrono::system_clock::now();
  std::cout << "FA:compute_stats" <<  (ms_t(end_time - start_time)).count() << "ms" << std::endl;

  std::cout << "ANI: " << stats.alignmentIdentity << std::endl;

  int seq_h_length = length(seq1);
  int seq_v_length = length(seq2);
  int seq_h_seed_length = (clippedEndPosition(row(align, 0)) - 1) -
                         clippedBeginPosition(row(align, 0));
  int seq_v_seed_length = (clippedEndPosition(row(align, 1)) - 1) -
                         clippedBeginPosition(row(align, 1));
  int seq_h_g_idx = 32816;
  int seq_v_g_idx = 50391;
  std::cout << stats.alignmentIdentity
               << "," << seq_h_length << "," << seq_v_length
               << "," << seq_h_seed_length << "," << seq_v_seed_length << std::endl;

  return 0;
}