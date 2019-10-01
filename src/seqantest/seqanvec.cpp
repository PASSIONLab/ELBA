// Created by Saliya Ekanayake on 2019-09-30.

#include <iostream>

#include <seqan/align.h>
#include <seqan/align_parallel.h>
#include <chrono>

typedef std::chrono::duration<double, std::milli> ms_t;

int main(int argc, char** argv){

  using TSequence = seqan::String<seqan::Dna>;
//  using TThreadModel = seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>;
  using TThreadModel = seqan::Serial;
  using TVectorSpec = seqan::Vectorial;
  using TExecPolicy = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;

  // dummy sequences
  std::string str1;
  std::string str2;
  str1.append(10000, 'A');
  str2.append(10000, 'A');

  seqan::StringSet<TSequence> seqs1;
  seqan::StringSet<TSequence> seqs2;

  appendValue(seqs1, TSequence{str1.c_str()});
  appendValue(seqs2, TSequence{str2.c_str()});

  TExecPolicy execPolicy;
  setNumThreads(execPolicy, 1);

  seqan::Score<int16_t, seqan::Simple> scoreAffine(1, -2, -1, -4);

  auto ts = std::chrono::system_clock::now();
  seqan::String<int16_t> scores = seqan::globalAlignmentScore(execPolicy, seqs1, seqs2, scoreAffine);
  auto te = std::chrono::system_clock::now();

  std::cout << "Score: " << scores[0] << "elapsed " << ((ms_t(te - ts)).count()) << "\n";

  return 0;
}

