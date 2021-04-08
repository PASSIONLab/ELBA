// Created by Saliya Ekanayake on 10/15/19 and modified by Giulia Guidi on 08/19/20.

#ifndef DIBELLA_CONTIGENTRY_HPP
#define DIBELLA_CONTIGENTRY_HPP

#include "Types.hpp"
#include "Defines.hpp"
#include "kmer/CommonKmers.hpp"

namespace dibella {

    void binary(ushort n, int* arr) 
    { 
        int nbit = 2;
        for(int i = 0; i < nbit; i++)
        { 
            arr[i] = n % 2; 
            n = n / 2; 
        }
    }

    struct ContigEntry {

        // GGGG: is dir/len consistent at this point? len can be derived from seq
        ushort dir;
        std::string seq;

        ContigEntry() : dir(4), seq("") {
        }

        // Overload + operator to add two ContigEntry objects
        ContigEntry operator+(const ContigEntry& other)
        {
            ContigEntry res;

            ushort rbit, lbit;
            ushort start, end;

            int mybin1[2] = {0, 0};
            int mybin2[2] = {0, 0};

            if(dir != 0) binary(dir, mybin1);
            if(other.dir != 0) binary(other.dir, mybin2);

            rbit = mybin1[0];
            lbit = mybin2[1];

            if(rbit != lbit)
            {
                start = mybin1[1]; 
                end   = mybin2[0]; 

                if(start == 0)
                {
                    if(end == 0) dir = 0;
                    else dir = 1;
                }
                else
                {
                    if(end == 0) dir = 2;
                    else dir = 3;      
                }
            }

            // Giulia keep an eye on this; esay buggy; these sequences need to be on the same strand when concatenating them
            if(dir == 1 || dir == 0) res.seq = seq + other.seq;
            else res.seq = other.seq + seq;

            res.dir = dir;
            return res;
        }

        // GGGG: I might not need this
        friend bool operator==(const ContigEntry& lhs, const ContigEntry& rhs)
        {
        	if(lhs.seq == rhs.seq) return true;
        	else return false;
        }
    };

    /* GGGG: matrix symmetrication removed */
    struct CeOutputHandler
    {
    	template <typename c, typename t>
    	void save(std::basic_ostream<c,t>  &os,
    			const dibella::ContigEntry &v,
    			uint64_t row,
    			uint64_t col)
    	{
    		int rc  = 0;
    		if(v.dir == 0 || v.dir == 3) rc = 1;

    		os << v.dir << "\t" << rc << "\t" << v.seq.length();
    	}
    };
}

#endif //DIBELLA_CONTIGENTRY_HPP
