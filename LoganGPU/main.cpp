#include <cmath>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <sstream>
#include <string>
#include <bits/stdc++.h>

#include "RunLoganAligner.hpp"

using namespace std;

void dnaSampleRun(string references, string queries, int xdrop, ushort kmerLen)
{
    vector<string> seqHs, seqVs;
    std::vector<SeedInterface> seeds;

    string   myInLine;
    ifstream rf(references);
    ifstream qf(queries);

    if(rf.is_open())
    {
        while(getline(rf, myInLine))
        {
            if(myInLine[0] == '>')
            {
                continue;
            }
            else
            {
                string seq = myInLine;
                seqHs.push_back(seq);

                SeedInterface seed(0, 0, kmerLen, kmerLen);
                seeds.push_back(seed);
            }
        }
        rf.close();
    }

    if(qf.is_open())
    {
        while(getline(qf, myInLine))
        {
            if(myInLine[0] == '>')
            {
                continue;
            }
            else
            {
                string seq = myInLine;
                seqVs.push_back(seq);
            }
        }
        qf.close();
    }

    std::vector<LoganResult> xscores(seqVs.size());
    RunLoganAlign(seqVs, seqHs, seeds, xscores, xdrop, kmerLen);

    std::cout << "Completed dnaSampleRun()" << std::endl;
}

int
main(int argc, char* argv[])
{
 	dnaSampleRun(argv[1], argv[2], stoi(argv[3]), stoi(argv[4]));
    return 0;
}
