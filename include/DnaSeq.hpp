#ifndef DNASEQ_H_
#define DNASEQ_H_

#include <cstdint>
#include <cstring>
#include <string>

class DnaSeq
{
public:
    static char    getcodechar(int c)  { return chartab[c]; }
    static uint8_t getcharcode(char c) { return codetab[(int)c]; }
    static char    getcharchar(char c) { return getcodechar(getcharcode(c)); }

    static constexpr char chartab[4+1] = {'A', 'C', 'G', 'T', '?'};

    static constexpr uint8_t codetab[256] =
    {
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
        4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };

    DnaSeq(char const *sequence, size_t len);
    DnaSeq(char const *sequence) : DnaSeq(sequence, strlen(sequence)) {}
    DnaSeq(std::string const& sequence) : DnaSeq(sequence.c_str(), sequence.size()) {}
    DnaSeq(const DnaSeq& rhs); //
    ~DnaSeq() { if (owns_memory) delete[] memory; }

    std::string ascii() const;
    size_t size() const { return 4*numbytes - remain; }
    size_t getnumbytes() const { return numbytes; }

    int operator[](size_t i) const; //
    bool operator<(const DnaSeq& rhs);
    bool operator==(const DnaSeq& rhs);
    bool operator!=(const DnaSeq& rhs) { return !(*this == rhs); }

    friend std::ostream& operator<<(std::ostream& stream, const DnaSeq& s)
    {
        stream << s.ascii();
        return stream;
    }

private:
    uint8_t *memory; /* each byte encodes up to 8 nucleotides */
    size_t numbytes; /* length of @memory */
    int remain; /* number of nucleotides in last byte */
    bool owns_memory; /* false if @memory pointer allocated somewhere else */
};

#endif
