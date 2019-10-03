// Created by Saliya Ekanayake on 2019-09-30.

#include <iostream>

#include <seqan/align.h>
#include <seqan/align_parallel.h>
#include <chrono>
#include <string>

typedef std::chrono::duration<double, std::milli> ms_t;

int main(int argc, char** argv){

  using TSequence = seqan::String<seqan::Dna>;
//  using TThreadModel = seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>;
  using TThreadModel = seqan::Serial;
  using TVectorSpec = seqan::Vectorial;
  using TExecPolicy = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;

  // dummy sequences

  seqan::Peptide seq1 = std::string(10000, 'A');
  seqan::Peptide seq2 = std::string(10000, 'A');

  seqan::StringSet<TSequence> seqs1;
  seqan::StringSet<TSequence> seqs2;

  TExecPolicy execPolicy;
  setNumThreads(execPolicy, 1);

  seqan::Score<int16_t, seqan::Simple> scoreAffine(1, -2, -1, -4);

  auto ts = std::chrono::system_clock::now();
  seqan::String<int16_t> scores = seqan::globalAlignmentScore(execPolicy, seq1, seq2, scoreAffine);
  auto te = std::chrono::system_clock::now();
  std::cout << "Vec Score: " << scores[0] << "elapsed " << ((ms_t(te - ts)).count()) << "\n";

  auto duration = ms_t(te- ts);


  seqan::Align<seqan::Peptide> align;
  resize(rows(align), 2);
  assignSource(row(align, 0), seq1);
  assignSource(row(align, 1), seq2);
  ts = std::chrono::system_clock::now();
  int score = seqan::globalAlignment(align, scoreAffine);
  te = std::chrono::system_clock::now();
  std::cout << "Just Score: " << scores[0] << "elapsed " << ((ms_t(te - ts)).count()) << "\n";

  duration += ms_t(te - ts);
  std::cout << "duration: " << duration.count() << std::endl;


  return 0;
}

