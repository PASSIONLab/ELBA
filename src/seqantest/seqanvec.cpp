// Created by Saliya Ekanayake on 2019-09-30.

#include <iostream>

#include <seqan/align.h>
#include <seqan/align_parallel.h>
#include <chrono>
#include <string>

#include <seqan/stream.h>  // for printint strings

typedef std::chrono::duration<double, std::milli> ms_t;

void vec1(){
  int N = 16;

  using TSequence = seqan::String<seqan::Dna>;
  using TThreadModel = seqan::WavefrontAlignment<seqan::BlockOffsetOptimization>;
//  using TThreadModel = seqan::Serial;
  using TVectorSpec = seqan::Vectorial;
  using TExecPolicy = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;
  seqan::Score<int16_t, seqan::Simple> scoreAffine(1, -2, -1, -4);

  // dummy sequences
//  seqan::Peptide seq1 = std::string(10000, 'A');
//  seqan::Peptide seq2 = std::string(10000, 'A');

  seqan::StringSet<TSequence> seqs1;
  seqan::StringSet<TSequence> seqs2;

  for (int i = 0; i < N; ++i) {
    std::string str1;
    std::string str2;
    str1.append(10000, 'A'+i);
    str2.append(10000, ('A'+i)-1);

    appendValue(seqs1, TSequence{str1.c_str()});
    appendValue(seqs2, TSequence{str2.c_str()});
  }


  TExecPolicy execPolicy;
  setNumThreads(execPolicy, 1);



  auto ts = std::chrono::system_clock::now();
  seqan::String<int16_t> scores = seqan::globalAlignmentScore(execPolicy, seqs1, seqs2, scoreAffine);
  auto te = std::chrono::system_clock::now();
  std::cout << "Vec Score: " << scores[0] << "elapsed " << ((ms_t(te - ts)).count()) << "\n";
  for (auto& s : scores){
    std::cout << s << " ";
  }
  std::cout << std::endl;

  auto duration = ms_t(te- ts);


  {
    ms_t duration;
    for (int i = 0; i < N; ++i) {
      std::string str1;
      std::string str2;
      str1.append(10000, 'A');
      str2.append(10000, ('A'+i)-1);

      seqan::Align<seqan::String<seqan::Dna>> align;
      resize(rows(align), 2);
      assignSource(row(align, 0), str1);
      assignSource(row(align, 1), str2);
      auto ts = std::chrono::system_clock::now();
      int score = seqan::globalAlignment(align, scoreAffine);
      auto te = std::chrono::system_clock::now();
      std::cout << "Just Score: " << score << "elapsed "
                << ((ms_t(te - ts)).count()) << "\n";

      duration += ms_t(te - ts);
    }
    std::cout << "duration: " << duration.count() << std::endl;
  }
}

void vec3(){
  using TSequence = seqan::String<seqan::AminoAcid>;
    using TThreadModel = seqan::Parallel;
    using TVectorSpec = seqan::Vectorial;
    using TExecPolicy = seqan::ExecutionPolicy<TThreadModel, TVectorSpec>;

    // dummy sequences
    TSequence seqH = "CGATT";
    TSequence seqV = "CGAAATT";

    seqan::StringSet<seqan::Gaps<TSequence>> seqs1;
    seqan::StringSet<seqan::Gaps<TSequence>> seqs2;

    for (size_t i = 0; i < 100; ++i)
    {
        appendValue(seqs1, seqan::Gaps<TSequence>(seqH));
        appendValue(seqs2, seqan::Gaps<TSequence>(seqV));
    }

    TExecPolicy execPolicy;
    setNumThreads(execPolicy, 4);

    seqan::Score<int16_t, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecBlosum62>> scoreAffine(-1, -4);

    seqan::String<int16_t> scores = seqan::globalAlignment(execPolicy, seqs1, seqs2, scoreAffine);

    for (int16_t score : scores)
        std::cout << "Score: " << score << "\n";

    std::cout << "SeqH: " << seqs1[0] << "\n";
    std::cout << "SeqV: " << seqs2[0] << "\n";
}


void vec2(){
  int N = 16;
  using TSequence = seqan::String<seqan::Dna>;
  seqan::StringSet<seqan::DnaString> collection1;
  seqan::StringSet<seqan::DnaString> collection2;

  seqan::StringSet<seqan::Align<seqan::DnaString>> ac;

  seqan::Score<int16_t, seqan::Simple> scoreAffine(1, -2, -1, -4);
  for (int i = 0; i < N; ++i) {
    std::string str1;
    std::string str2;
    str1.append(10000, 'A'+i);
    str2.append(10000, ('A'+i)-1);
    seqan::Align<seqan::Peptide> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), str1);
    assignSource(row(align, 1), str2);

//    appendValue(collection1, TSequence{str1.c_str()});
//    appendValue(collection2, TSequence{str2.c_str()});
    appendValue(ac, align);
  }


//  assert(seqan::length(collection1) == seqan::length(collection2));

  auto ts = std::chrono::system_clock::now();
//  seqan::String<int> scores =  seqan::globalAlignmentScore(collection1, collection2, scoreAffine);
  seqan::String<int> scores =  seqan::globalAlignment(ac, scoreAffine);
  auto te = std::chrono::system_clock::now();
  std::cout << "Vec Score: " << scores[0] << "elapsed " << ((ms_t(te - ts)).count()) << "\n";

  assert(seqan::length(scores) == seqan::length(collection1));
  
  /*{
    ms_t duration;
    for (int i = 0; i < N; ++i) {
      std::string str1;
      std::string str2;
      str1.append(10000, 'A');
      str2.append(10000, ('A'+i)-1);

      TSequence seq1{str1.c_str()};
      TSequence seq2{str2.c_str()};
//      seqan::Align<seqan::String<seqan::Dna>> align;
//      resize(rows(align), 2);
//      assignSource(row(align, 0), str1);
//      assignSource(row(align, 1), str2);
      auto ts = std::chrono::system_clock::now();
      int score = seqan::globalAlignmentScore(seq1, seq2, scoreAffine);
      auto te = std::chrono::system_clock::now();
      std::cout << "Just Score: " << score << "elapsed "
                << ((ms_t(te - ts)).count()) << "\n";

      duration += ms_t(te - ts);
    }
    std::cout << "duration: " << duration.count() << std::endl;
  }*/

}

int main(int argc, char** argv){
//  vec2();
  vec3();
  return 0;
}

