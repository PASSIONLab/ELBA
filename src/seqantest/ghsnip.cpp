#include <iostream>

#include <seqan/align_parallel.h>
#include <seqan/stream.h>  // for printint strings

int main()
{
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

    return EXIT_SUCCESS;
}
