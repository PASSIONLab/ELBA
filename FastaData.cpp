//
// Created by Saliya Ekanayake on 1/7/19.
//

#include <iostream>
#include <mpi.h>
#include "FastaData.hpp"

FastaData::~FastaData() {

}

FastaData::FastaData(std::shared_ptr<char> data, int l_start, int l_end) {
  this->data = data;
  this->l_start = l_start;
  this->l_end = l_end;

  id_starts = new std::vector<int>();
  seq_starts = new std::vector<int>();

  l_seq_count = 0;
  char c;
  bool in_name = false;
  /* Assume the FASTA content is valid */
  for (int i = l_start; i <= l_end; ++i) {
    c = data.get()[i];
    if (c == '>') {
      id_starts->push_back(i);
      in_name = true;
      ++l_seq_count;
    }
    if (c == '\n' && in_name && i + 1 <= l_end) {
      seq_starts->push_back(i + 1);
      in_name = false;
    }
  }

  g_seq_offset = 0;
  MPI_Exscan(&l_seq_count, &g_seq_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

}

void FastaData::print() {
  // Test print
  for (int i = l_start; i <= l_end; ++i) {
    std::cout << data.get()[i];
  }

  std::cout<<"\n===========\n";
  std::cout<<"l_seq_count: "<<l_seq_count<<" g_seq_offset: "<<g_seq_offset<<"\n";
  for (int i = 0; i < id_starts->size(); ++i){
    char *beg = data.get() + (*id_starts)[i];
    char *end = data.get() + ((*seq_starts)[i] - 1);
    std::cout.write(beg, end - beg);
    std::cout<<"\n"<<"(end - beg)="<<(end-beg)<<" end_char=|"<<data.get()[((*seq_starts)[i] - 2)]<<"|\n";
    beg = data.get() + (*seq_starts)[i];
    end = data.get() + (i+1 < id_starts->size() ? ((*id_starts)[i+1] - 1) : l_end+1);
    std::cout.write(beg, end - beg);
    std::cout<<"\n";
  }
  std::cout<<"\n";
}

