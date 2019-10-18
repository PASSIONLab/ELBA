// Created by Saliya Ekanayake on 10/18/19.
#include <iostream>
#include "include/kmer/KmerOps.hpp"


std::vector<std::string> func(
    std::string& w1, std::string&w2, int k, std::string alph,
    std::string& prefix, int level, bool recurring){
  bool limit = prefix[prefix.length() - 1] == w2[level - 1];
  if (level == k){
    std::vector<std::string> wlist;
    if (recurring){
      for (auto& let : alph){
        wlist.push_back(prefix+let);
      }
    } else {
      if(prefix == w2.substr(0, level - 1)){
//        for (char let = alph[w1[level - 1]])
      }
    }
  }
}

void genStrings(){
  int k = 4;
  std::string alph("ABCD");

  std::string w1("BCA");
  std::string w2("CCB");


}

/*! NOOP, this is not the same that I want */
std::vector<std::string> generateAllPossibleStrings(std::string start, std::string end, int k) {

  std::vector<std::string> variants;
//  char startArray[] = start.toCharArray();
//  char endArray[] = end.toCharArray();
//  char currentArray[] = Arrays.copyOf(startArray, startArray.length);
  std::string current_str = start;
  variants.push_back(start);

  //We check if the start std::string is really above the end std::string as specified
  //We output an empty std::string if it is not the case
  bool possible = true;
  for(int i = 0; i<k; i++)
    possible = possible && (start[i]<=end[i]);
  if (!possible)
    return variants;


  while(end != current_str){
    current_str[k-1]+=1;
    int i = k-1;
    while(current_str[i]>end[i]){
      current_str[i]=start[i];
      i--;
      current_str[i]++;
    }
    variants.push_back(current_str);
  }

  return variants;
}

int main(int argc, char** argv){
//  auto parops = ParallelOps::init(&argc, &argv);
//  ushort k = 2;
//  Alphabet alph(Alphabet::PROTEIN);
//  pisa::KmerOps::generate_S(k, alph, parops);
//  parops->eardown_parallelism();

//  std::string start("MAC");
//  std::string end("PQN");
//  int k = 3;
//  auto nbrs = generateAllPossibleStrings(start, end, k);
//  for (auto& str : nbrs){
//    std::cout << str << std::endl;
//  }

  const char* seq = "ABCDEF";
  std::string s(seq+1, 3);
  std::cout << s << std::endl;
}
