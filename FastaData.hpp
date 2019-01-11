//
// Created by Saliya Ekanayake on 1/7/19.
//

#ifndef LBL_DAL_FASTADATA_HPP
#define LBL_DAL_FASTADATA_HPP


#include <cstdlib>
#include <memory>
#include <vector>

class FastaData {
public:
  FastaData(std::shared_ptr<char> data, int l_start, int l_end);
  ~FastaData();

  void print();

private:
  std::shared_ptr<char> data;
  int l_seq_count;
  int g_seq_offset;
  int l_start, l_end;
  std::vector<int>* id_starts = nullptr;
  std::vector<int>* seq_starts = nullptr;


};


#endif //LBL_DAL_FASTADATA_HPP
