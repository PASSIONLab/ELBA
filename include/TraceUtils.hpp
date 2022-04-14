// Created by Saliya Ekanayake on 2/10/19.

#ifndef LBL_DAL_DEBUGUTILS_HPP
#define LBL_DAL_DEBUGUTILS_HPP

#include "ParallelOps.hpp"
#include "Types.hpp"

struct TimePod {
  std::unordered_map<std::string, ticks_t> times;
  std::string names[21] = {"Main",
                           "Main:newDFD()",
                           "Dfd:PfrReadFasta()",
                           "Dfd:newFD()",
                           "KmerOp:GenerateA:CardinalityHLL()",
                           "KmerOp:GenerateA:FirstPass()",
                           "KmerOp:GenerateA:SecondPass()",
                           "KmerOp:GenerateA:SpMatA()",
                           "Main:GenerateA()",
                           "Main:At()",
                           "Main:AAt()",
                           "Main:DfdWait()",
                           "Dfd:MPI_Waitall(seqs)",
                           "Main:DprAlign()",
                           "Main:TransitiveReduction()",
                           "Main:ExtractContig()",
                           "CreateContig:GetRead2Contigs()",
                           "CreateContig:GetRead2ProcAssignments()",
                           "CreateContig:InducedSubgraphs2Procs()",
                           "CreateContig:ReadExchange()",
                           "CreateContig:LocalAssembly()"
  };

  std::string to_string() {
    std::string str = "\nINFO: Program timings:\n";
    ticks_t ts, te;
    for (const auto &name : names) {
      if (times.find("Start" + name) != times.end()) {
        ts = times["Start" + name];
        te = times["End" + name];
        str.append("  ").append(name).append(":")
          .append(std::to_string((ms_t(te - ts)).count())).append(" ms\n");
      } else {
        str.append("  ").append(name).append(" not evaluated.\n");
      }
    }
    return str;
  }
};

class TraceUtils {
public:
  //  static void print_fasta_data(const std::unique_ptr<FastaData> &Dfd,
  //  const std::shared_ptr<ParallelOps> &parops);
  static void print_msg(const std::string &title, const std::string &msg,
                        const std::shared_ptr<ParallelOps> &parops);

  static void
  print_msg_on_rank(const std::string &title, const std::string &msg,
                    const std::shared_ptr<ParallelOps> &parops,
                    int rank);

  explicit TraceUtils(bool is_print_rank);

  void print_str(std::string str);

private:
  bool is_print_rank;
};

//void _debug_print_fasta_data(const std::unique_ptr<FastaData> &fd,
//                             const std::shared_ptr<ParallelOps> &p_ops){
//  /* Debug print of FastData instance*/
//  int flag;
//  if (p_ops->world_proc_rank > 0)
//    MPI_Recv(&flag, 1, MPI_INT, p_ops->world_proc_rank-1,
//             99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//  std::cout<<"\nRank: "<<p_ops->world_proc_rank<<"\n------------------------\n";
//  fd->print();
//  if (p_ops->world_proc_rank < p_ops->world_procs_count - 1) {
//    MPI_Send(&flag, 1, MPI_INT, p_ops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//}

//void _print_str(const std::shared_ptr<ParallelOps> &p_ops, const std::string &str){
//  int flag;
//  if (p_ops->world_proc_rank > 0)
//    MPI_Recv(&flag, 1, MPI_INT, p_ops->world_proc_rank-1,
//             99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//  std::cout<<"\nRank: "<<p_ops->world_proc_rank<<"\n------------------------\n";
//  std::cout<<str;
//  if (p_ops->world_proc_rank < p_ops->world_procs_count - 1) {
//    MPI_Send(&flag, 1, MPI_INT, p_ops->world_proc_rank + 1, 99, MPI_COMM_WORLD);
//  }
//  MPI_Barrier(MPI_COMM_WORLD);
//}

#endif //LBL_DAL_DEBUGUTILS_HPP
