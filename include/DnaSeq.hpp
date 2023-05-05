#ifndef DNASEQ_H_
#define DNASEQ_H_

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>

/*
 * Originally I wanted this class to be able to be used under two different
 * circumstances: (1) where a given DnaSeq object is responsible for allocating
 * its buffer memory and freeing it and (2) where a DnaBuffer object provides
 * the memory and the DnaSeq object provided a pointer into that buffer. In
 * C this would be easy to do, but I don't fully understand C++ yet and so
 * ran into problems where the copy constructor was being automatically called
 * in certain cases. The copy constructor was implemented to allocate a new
 * memory buffer. This in itself is not problematic, but when I wanted to use
 * the raw pointer of a given DnaSeq object to find its offset within a DnaBuffer,
 * I was getting garbage because the the DnaSeq pointer was being copied.
 *
 * Because I don't want to spend any more time figuring this out, and more importantly,
 * because I actually don't need to have two different storage paradigms for
 * ELBA, I'm just going to get rid of circumstance (1) listed above.
 *
 * Therefore, this DnaSeq object is currently to be used ONLY in a context
 * where its memory is given by an underlying DnaBuffer object. I may change
 * this later but right now I don't see any reason to.
 */

class DnaSeq
{
public:
    /*
     * Null constructor for resizing vectors.
     */
    DnaSeq() : len(0), memory(nullptr) {}

    DnaSeq(size_t len, uint8_t *mem) : len(len), memory(mem) {}

    /*
     * char const *s - ASCII DNA sequence. Can include lower case letters. Ns converted to As.
     *                 Non-nucleotide characters cause undefined behaviour.
     * size_t len    - sequence length (number of nucleotides)
     * uint8_t* mem  - pointer to where compressed sequence should be written to and stored.
     */
    DnaSeq(char const *s, size_t len, uint8_t *mem) : DnaSeq(len, mem) { compress(s); }

    /*
     * Copy constructor. Copied DnaSeq objects reference
     * the same memory buffer.
     */
    DnaSeq(const DnaSeq& rhs) : len(rhs.len), memory(rhs.memory) {}

    /*
     * Get an ASCII copy of the sequence (upper-case).
     */
    std::string ascii() const;

    /*
     * Get number of nucleotides in sequence.
     */
    size_t size() const { return len; }

    /*
     * Get the number of bytes used to encode the sequence.
     */
    size_t numbytes() const { return bytesneeded(size()); }

    /*
     * Because we can encode up to 4 nucleotides in a single byte,
     * we will have some leftover unused bits in the last byte of
     * the sequence if the number of nucleotides is not divisible
     * by 4. This returns the number of nucleotides that aren't
     * encoded in the last byte. will be a value in the range [0..3].
     */
    int remainder() const { return 4*numbytes() - size(); }

    /*
     * Returns raw pointer to the first byte of memory.
     */
    const uint8_t* data() const { return memory; }

    /*
     * Copies data into location (meant to be used with DnaBufferContainer).
     */
    void copyto(size_t *readlen, uint8_t *mem) const
    {
        *readlen = len;
        std::memcpy(mem, memory, numbytes());
    }

    /*
     * returns integer code of given position:
     * 0 means A, 1 means C, 2 means G, and 3 means T.
     */
    int operator[](size_t i) const;

    /*
     * Lexicographical comparison operator.
     */
    bool operator<(const DnaSeq& rhs);

    /*
     * Equality comparision operators.
     */
    bool operator==(const DnaSeq& rhs);
    bool operator!=(const DnaSeq& rhs) { return !(*this == rhs); }

    /*
     * same as operator[].
     */
    int regular_at(size_t i) const { return (*this)[i]; }

    /*
     * same as operator[], except it interprets the sequence
     * as its reverse complement.
     */
    int revcomp_at(size_t i) const { return 3 - (*this)[len-1-i]; }

    /*
     * Returns the number of bytes needed to encode a given number of nucleotides.
     */
    static size_t bytesneeded(size_t n) { return (n+3)/4; }

    /*
     * getcodechar : [0,1,2,3,4] -> [A,C,G,T,X]
     * getcharcode : [A,a,C,c,G,g,T,t,N,n,...] -> [0,0,1,1,2,2,3,3,0,0,X]
     * getcharchar : [A,a,C,c,G,g,T,t,N,n,...] -> [A,A,C,C,G,G,T,T,A,A,X]
     */
    static char    getcodechar(int c)  { return chartab[c]; }
    static uint8_t getcharcode(char c) { return codetab[(int)c]; }
    static char    getcharchar(char c) { return getcodechar(getcharcode(c)); }

    static constexpr char chartab[4+1] = {'A', 'C', 'G', 'T', 'X'};
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

    friend std::ostream& operator<<(std::ostream& stream, const DnaSeq& s)
    {
        stream << s.ascii();
        return stream;
    }

private:
    size_t len; /* number of nucleotides encoded starting at @memory */
    uint8_t *memory; /* each byte encodes up to 4 nucleotides */

    /*
     * This is to be used only by the constructors of DnaSeq!
     */
    void compress(char const *s);
};

#endif
