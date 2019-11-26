// Created by Saliya Ekanayake on 11/26/19.
#include <iostream>
#include <string>
#include <seqan/align.h>

typedef std::chrono::duration<double, std::milli> ms_t;

int main(int argc, char** argv){
  std::string str_255("ssqlipnispdsftvaastgmlsgkshemlydaetgrkisqldwkiknvailkgdiswdpysfltlnargwtslasgsgnmddyawmnenqsewtdhsshpatnvnhaneydlnvkgwllqdenykagitagyqetrfswtatggsysynngaytgnfpkgvrvigynqrfsmpyiglagqyrindfelnalfkfsdwvrahdndehymrdltfrektsgsryygtvinagyyvtpnakvfaeftyskydegkggtqtidknsgdsvsiggdaagisnknytvtaglqyrf");
  std::string str_288("pyieifeqprqrgmrfrykcegrsagsipgehstdnnktfpsiqilnyfgkvkirttlvtknepykphphdlvgkdcrdgyyeaefgperrvlsfqnlgiqcvkkkdlkesislriskkinpfnvpeeqlhnideydlnvvrlcfqaflpdehgnytlalpplisnpiydnrapn");

  seqan::Peptide seq_255(str_255);
  seqan::Peptide seq_288(str_288);

  int gap_ext = -2;
  int gap_open = -11;
  seqan::Blosum62 blosum62(gap_ext, gap_open);

  seqan::Align<seqan::Peptide> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq_255);
  assignSource(row(align, 1), seq_255);

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
  
  return 0;
}