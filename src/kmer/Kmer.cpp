#include <bitset>
#include <string>
#include <iostream>

#include "../../include/kmer/Kmer.hpp"

#ifdef __APPLE__
#include <machine/endian.h>
#define H2BE htonll
#define BE2H ntohll
#else
#include <endian.h>
#define H2BE htobe64
#define BE2H be64toh
#endif

using namespace std;

static const uint64_t twin_table[256] = {
0xFF, 0xBF, 0x7F, 0x3F, 0xEF, 0xAF, 0x6F, 0x2F,
0xDF, 0x9F, 0x5F, 0x1F, 0xCF, 0x8F, 0x4F, 0x0F,
0xFB, 0xBB, 0x7B, 0x3B, 0xEB, 0xAB, 0x6B, 0x2B,
0xDB, 0x9B, 0x5B, 0x1B, 0xCB, 0x8B, 0x4B, 0x0B,
0xF7, 0xB7, 0x77, 0x37, 0xE7, 0xA7, 0x67, 0x27,
0xD7, 0x97, 0x57, 0x17, 0xC7, 0x87, 0x47, 0x07,
0xF3, 0xB3, 0x73, 0x33, 0xE3, 0xA3, 0x63, 0x23,
0xD3, 0x93, 0x53, 0x13, 0xC3, 0x83, 0x43, 0x03,
0xFE, 0xBE, 0x7E, 0x3E, 0xEE, 0xAE, 0x6E, 0x2E,
0xDE, 0x9E, 0x5E, 0x1E, 0xCE, 0x8E, 0x4E, 0x0E,
0xFA, 0xBA, 0x7A, 0x3A, 0xEA, 0xAA, 0x6A, 0x2A,
0xDA, 0x9A, 0x5A, 0x1A, 0xCA, 0x8A, 0x4A, 0x0A,
0xF6, 0xB6, 0x76, 0x36, 0xE6, 0xA6, 0x66, 0x26,
0xD6, 0x96, 0x56, 0x16, 0xC6, 0x86, 0x46, 0x06,
0xF2, 0xB2, 0x72, 0x32, 0xE2, 0xA2, 0x62, 0x22,
0xD2, 0x92, 0x52, 0x12, 0xC2, 0x82, 0x42, 0x02,
0xFD, 0xBD, 0x7D, 0x3D, 0xED, 0xAD, 0x6D, 0x2D,
0xDD, 0x9D, 0x5D, 0x1D, 0xCD, 0x8D, 0x4D, 0x0D,
0xF9, 0xB9, 0x79, 0x39, 0xE9, 0xA9, 0x69, 0x29,
0xD9, 0x99, 0x59, 0x19, 0xC9, 0x89, 0x49, 0x09,
0xF5, 0xB5, 0x75, 0x35, 0xE5, 0xA5, 0x65, 0x25,
0xD5, 0x95, 0x55, 0x15, 0xC5, 0x85, 0x45, 0x05,
0xF1, 0xB1, 0x71, 0x31, 0xE1, 0xA1, 0x61, 0x21,
0xD1, 0x91, 0x51, 0x11, 0xC1, 0x81, 0x41, 0x01,
0xFC, 0xBC, 0x7C, 0x3C, 0xEC, 0xAC, 0x6C, 0x2C,
0xDC, 0x9C, 0x5C, 0x1C, 0xCC, 0x8C, 0x4C, 0x0C,
0xF8, 0xB8, 0x78, 0x38, 0xE8, 0xA8, 0x68, 0x28,
0xD8, 0x98, 0x58, 0x18, 0xC8, 0x88, 0x48, 0x08,
0xF4, 0xB4, 0x74, 0x34, 0xE4, 0xA4, 0x64, 0x24,
0xD4, 0x94, 0x54, 0x14, 0xC4, 0x84, 0x44, 0x04,
0xF0, 0xB0, 0x70, 0x30, 0xE0, 0xA0, 0x60, 0x20,
0xD0, 0x90, 0x50, 0x10, 0xC0, 0x80, 0x40, 0x00
};

// use:  km = Kmer();
// pre:  
// post: the DNA string in km is AA....AAA (k times A) 
Kmer::Kmer()
{
    for (size_t i = 0; i < N_LONGS; i++)
    {
        longs[i] = 0;
    }
}

// use:  _km = Kmer(km);
// pre:  s[0],...,s[k] are all equal to 'A','C','G' or 'T'
// post: the DNA string in _km and is the same as in km 
Kmer::Kmer(const Kmer& o)
{
    for (size_t i = 0; i < N_LONGS; i++)
    {
        longs[i] = o.longs[i];
    }
}

// use:  km = Kmer(s);
// pre:  s[0],...,s[k] are all equal to 'A','C','G' or 'T'
// post: the DNA string in km is now the same as s
Kmer::Kmer(const char *s) {
    set_kmer(s);
}

// use:  _km = km;
// pre:  
// post: the DNA string in _km and is the same as in km 
Kmer& Kmer::operator=(const Kmer& o)
{
    if (this != &o)
    {
        for (size_t i = 0; i < N_LONGS; i++)
        {
          longs[i] = o.longs[i];
        }
    }
    return *this;
}

// use:  km = Kmer();
// pre:  
// post: The last bit in the bit array which stores the DNA string has been set to 1
//       which indicates that the km is invalid
void Kmer::set_deleted()
{
    memset(bytes.data(),0xff,N_BYTES);
}

// use:  b = (km1 < km2);
// pre:  
// post: b is true <==> the DNA strings in km1 is alphabetically smaller than
//                      the DNA string in km2 
bool Kmer::operator<(const Kmer& o) const
{
    bool r = false;
    for (size_t i = 0; i < N_LONGS; ++i)
    {
        if (longs[i] < o.longs[i])
            return true;
        if (longs[i] > o.longs[i])
            return false;
    }
    return false;
}

// use:  b = (km1 == km2);
// pre:  
// post: b is true <==> the DNA strings in km1 and km2 are equal
bool Kmer::operator==(const Kmer& o) const
{
    for (size_t i = 0; i < N_LONGS; i++)
    {
        if (longs[i] != o.longs[i])
        {
           return false;
        }
    }
    return true;
}

// use:  km.set_kmer(s);
// pre:  s[0],...,s[k-1] are all 'A','C','G' or 'T'
// post: The DNA string in km is now equal to s 
void Kmer::set_kmer(const char *s)  {
    size_t i,j,l;
    memset(bytes.data(),0,N_BYTES);
        
    for (i = 0; i < k; ++i)
    {
        j = i % 32;
        l = i/32;
        assert(*s != '\0');
        
        size_t x = ((*s) & 4) >> 1;
        longs[l] |= ((x + ((x ^ (*s & 2)) >>1)) << (2*(31-j)));    
        s++; 
    }
}

// This returns a vector of kmers from a long sequence
// This does not check for Ns or choose the lesser representation
std::vector<Kmer> Kmer::getKmers(std::string seq)
{
    for (auto & c: seq) c = toupper(c); // ensure uppercase
    if (seq.size() < k) return std::vector<Kmer>();

    int bufsize = std::max((int) N_LONGS, (int) (seq.size()+31)/32) + 2;
    int numLongs = (k+31)/32;
    assert(numLongs <= N_LONGS);
    int lastLong = numLongs - 1;
    assert(lastLong >= 0 && lastLong < N_LONGS);

    std::vector<Kmer> kmers(seq.size() - k + 1, Kmer());
    uint64_t buf[ bufsize ];
    uint8_t *bufPtr = (uint8_t*) buf;
    memset(buf, 0, bufsize*8);
    const char *s = seq.c_str();

    // calculate binary along entire sequence
    for (int i = 0; i < seq.size(); ++i)
    {
        int j = i % 32;
        int l = i/32;
        assert(*s != '\0');
        
        size_t x = ((*s) & 4) >> 1;
        buf[l] |= ((x + ((x ^ (*s & 2)) >>1)) << (2*(31-j)));
        s++;
    }

    // fix to big endian
    for (int l = 0; l < bufsize; l++)
    {
        buf[l] = H2BE(buf[l]);
    }
    const uint64_t mask = ((int64_t) 0x3);
    uint64_t endmask = 0;

    if (k%32)
    {
        endmask = (((uint64_t) 2) << (2*(31-(k%32)) + 1)) - 1;

    }
    endmask = ~endmask;

    // make 4 passes for each phase of 2 bits per base
    for (int shift = 0; shift < 4; shift++)
    {
        if (shift > 0)
        {
            // shift all bits by 2 in the buffer of longs
            for(int l = 0; l < bufsize - 1; l++)
            {
                buf[l] = (~mask & (BE2H(buf[l])<<2)) | (mask & (BE2H(buf[l+1]) >> 62));
                buf[l] = H2BE(buf[l]);
            }
        }

        // enumerate the kmers in the phase
        for(int i = shift; i < kmers.size(); i += 4)
        {
            int byteOffset = i/4;
            assert(byteOffset+numLongs*8 <= bufsize*8);

            for(int l = 0; l < numLongs; l++)
            {
                    kmers[i].longs[l] = BE2H( *( (uint64_t *) (bufPtr + byteOffset + l*8) ) );
            }

            // set remaining bits to 0
            kmers[i].longs[lastLong] &= endmask;
            
            for(int l = numLongs; l < N_LONGS; l++)
            {
                    kmers[i].longs[l] = 0;
            }

    #ifdef DBUG
            Kmer test(seq.c_str() + i);
            if(kmers[i] != test)
            {
                fprintf(stderr, "Invalid kmer kmer %d %s %s from %s\n", i, kmers[i].toString().c_str(), test.toString().c_str(), seq.c_str());
            }
    #endif
        }
    }
    return kmers; 
}

// use:  i = km.hash();
// pre:   
// post: i is the hash value of km 
uint64_t Kmer::hash() const {
    uint64_t ret;
    assert(k_bytes>0);
    return MurmurHash3_x64_64((const void*)bytes.data(), N_BYTES);
}

// use:  rep = km.rep();
// pre:   
// post: rep is km.twin() if the DNA string in km.twin() is alphabetically smaller than
//       the DNA string in km, else rep is km
Kmer Kmer::rep() const {
  Kmer tw = twin();
  return (tw < *this) ? tw : *this;
}

// use:  tw = km.twin();
// pre:   
// post: tw is the twin kmer with respect to km,
//       i.e. if the DNA string in km is 'GTCA'
//          then the DNA string in tw is 'TGAC'
Kmer Kmer::twin() const {
  
    Kmer km(*this);
    size_t nlongs = (k+31)/32;
    
    for (size_t i = 0; i < nlongs; i++)
    {
        uint64_t v = longs[i];
        km.longs[nlongs-1-i] =  
          (twin_table[v & 0xFF] << 56) | 
          (twin_table[(v>>8) & 0xFF] << 48) | 
          (twin_table[(v>>16) & 0xFF] << 40) | 
          (twin_table[(v>>24) & 0xFF] << 32) |
          (twin_table[(v>>32) & 0xFF] << 24) | 
          (twin_table[(v>>40) & 0xFF] << 16) | 
          (twin_table[(v>>48) & 0xFF] << 8)  |
          (twin_table[(v>>56)]);
    }

    size_t shift = (k%32) ? 2*(32-(k%32)) : 0 ;
    uint64_t shiftmask = (k%32) ? (((1ULL<< shift)-1) << (64-shift)) : 0ULL;
  
    km.longs[0] = km.longs[0] << shift;
    for (size_t i = 1; i < nlongs; i++)
    {
        km.longs[i-1] |= (km.longs[i] & shiftmask) >> (64-shift);
        km.longs[i] = km.longs[i] << shift;    
    }
  
  return km;
}

// use:  link = km.getLink(index);
// pre:  0 <= index < 8
// post: gives the forward kmer with the (index % 4) character in 'A','C','G' or 'T' if index < 4
//       else the backward kmer with the (index % 4) character in 'A','C','G' or 'T'
Kmer Kmer::getLink(const size_t index) const
{
    assert(index >= 0 && index < 8);
    char c;

    switch (index % 4)
    {
        case 0: c = 'A'; break;
        case 1: c = 'C'; break;
        case 2: c = 'G'; break;
        case 3: c = 'T'; break;
    }
    
    return (index < 4) ? forwardBase(c) : backwardBase(c);
}

// use:  fw = km.forwardBase(c)
// pre:  
// post: fw is the forward kmer from km with last character c,
//       i.e. if the DNA string in km is 'ACGT' and c equals 'T' then
//       the DNA string in fw is 'CGTT'
Kmer Kmer::forwardBase(const char b) const
{
    Kmer km(*this);

    km.longs[0] = km.longs[0] << 2;
    size_t nlongs = (k+31)/32;

    for (size_t i = 1; i < nlongs; i++)
    {
        km.longs[i-1] |= (km.longs[i] & (3ULL<<62)) >> 62;
        km.longs[i]  = km.longs[i] << 2;
    }

    uint64_t x = (b & 4) >>1;
    km.longs[nlongs-1] |= (x + ((x ^ (b & 2)) >>1 )) << (2*(31-((k-1)%32)));

    return km;
}


bool Kmer::equalUpToLastBase(const Kmer & rhs)
{
	Kmer Ashifted = rhs.backwardBase('A');
	Kmer Bshifted = backwardBase('A');
	return (Ashifted == Bshifted);
}

// use:  bw = km.backwardBase(c)
// pre:  
// post: bw is the backward kmer from km with first character c,
//       i.e. if the DNA string in km is 'ACGT' and c equals 'T' then
//       the DNA string in bw is 'TACG'
Kmer Kmer::backwardBase(const char b) const
{
    Kmer km(*this);
    
    size_t nlongs = (k+31)/32;
    km.longs[nlongs-1] = km.longs[nlongs-1] >>2;
    km.longs[nlongs-1] &= (k%32) ? (((1ULL << (2*(k%32)))-1) << 2*(32-(k%32))) : ~0ULL;

    for (size_t i = 1; i < nlongs; i++)
    {
        km.longs[nlongs-i] |= (km.longs[nlongs-i-1] & 3ULL) << 62;
        km.longs[nlongs-i-1] = km.longs[nlongs-i-1] >>2;
    }
    uint64_t x = (b & 4) >> 1;
    km.longs[0] |= (x + ((x ^ (b & 2)) >> 1)) << 62;

    return km;
}

// use:  km.printBinary();
// pre:   
// post: The bits in the binary representation of the 
//       DNA string for km has been printed to stdout
std::string Kmer::getBinary() const
{
    size_t nlongs = N_LONGS;
    std::string r;
    r.reserve(64*nlongs);
    for (size_t i = 0; i < nlongs; i++)
    {
        r.append(std::bitset<64>(longs[i]).to_string<char,std::char_traits<char>,std::allocator<char>>());
    }
    return r;
}

// use:  km.toString(s);
// pre:  s has space for k+1 elements
// post: s[0,...,k-1] is the DNA string for the Kmer km and s[k] = '\0'
void Kmer::toString(char * s) const
{
    size_t i, j, l;
    for (i = 0; i < k; i++)
    {
        j = i % 32;
        l = i / 32;

        switch(((longs[l]) >> (2*(31-j)) )& 0x03 )
        {
            case 0x00: *s = 'A'; ++s; break;
            case 0x01: *s = 'C'; ++s; break;
            case 0x02: *s = 'G'; ++s; break;
            case 0x03: *s = 'T'; ++s; break;  
        }
    }
    *s = '\0';
}

std::string Kmer::toString() const
{
    char buf[MAX_K];
    toString(buf);
    return std::string(buf);
}

// use:  set_k(k);
// pre:  this method has not been called before and 0 < k < MAX_K
// post: The Kmer size has been set to k
void Kmer::set_k(unsigned int _k)
{
    assert(_k < MAX_K);
    assert(_k > 0);
    if(_k == k)
    {
        return; // ok to call more than once
    }

    k = _k;
    k_bytes = (_k+3)/4;
    k_modmask = (1 << (2*((k%4)?k%4:4)) )-1;
}

unsigned int Kmer::k = 0;
unsigned int Kmer::k_bytes = 0;
unsigned int Kmer::k_modmask = 0;