//
// Created by Saliya Ekanayake on 1/6/19.
//

#include <iostream>
#include <fstream>
#include <memory>
#include <cassert>
#include <limits>
#include "ParallelFastaReader.hpp"
#include "FastaData.hpp"

void ParallelFastaReader::readFasta(const char *file, ushort overlap,
                                    int rank, int world_size,
                                    std::shared_ptr<FastaData> &fd) {
  MPI_File f;

  int err = MPI_File_open(MPI_COMM_WORLD, file,
      MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
  if (err) {
    if (rank == 0) fprintf(stderr, "Couldn't open file %s\n", file);
    MPI_Finalize();
    exit(2);
  }

  /* The idea:
   * if sequence names and contents are roughly equal
   * then dividing the file size by the number of ranks
   * should land close to l_seq_count sequences. We'll fix
   * any mismatch later
   */

  /* Thanks, Rob Lathem for making my life a bit easy here
   * https://stackoverflow.com/a/12942718 */
  MPI_Offset g_start;
  int64_t l_chunk_size;

  /* read in relevant chunk of file into "chunk",
   * which starts at location in the file g_start
   * and has size l_chunk_size
   */

  MPI_Offset g_end;
  MPI_Offset file_size;

  /* figure out who reads what */
  MPI_File_get_size(f, &file_size);
  file_size--;  /* get rid of text file eof */
  l_chunk_size = static_cast<int64_t>(file_size / world_size);
  g_start = rank * l_chunk_size;
  g_end = g_start + l_chunk_size - 1;
  if (rank == world_size - 1) g_end = file_size - 1;

  /* add overlap to the end of everyone's chunk except last proc... */
  if (rank != world_size - 1)
    g_end += overlap;

  l_chunk_size = static_cast<int64_t>(g_end - g_start + 1);

  /* create shared pointer to data */
  std::shared_ptr<char> data(
      static_cast<char*>(malloc((l_chunk_size+1) * sizeof(char))), free);

  /* everyone reads in their part */
  char* chunk = data.get();
  // TODO - Saliya: fix if l_chunk_size > int max. Read multiple times
  assert(l_chunk_size <= std::numeric_limits<int>::max());
  MPI_File_read_at_all(f, g_start, chunk, static_cast<int>(l_chunk_size),
                       MPI_CHAR, MPI_STATUS_IGNORE);
  chunk[l_chunk_size] = '\0';


  /*
   * everyone calculate what their start and end *really* are by going
   * from the first newline after start to the first newline after the
   * overlap region starts (eg, after end - overlap + 1)
   */

  int64_t l_start = 0, l_end = l_chunk_size - 1;
  if (rank != 0) {
    while (chunk[l_start] != '>') l_start++;
  }
  if (rank != world_size - 1) {
    l_end -= overlap;
    while (chunk[l_end] != '>') l_end++;
    // minus 2 because we don't need '>' as well as the '\n' before that
    l_end -= 2;
  }

  fd = std::make_shared<FastaData>(data, l_start, l_end);
}



