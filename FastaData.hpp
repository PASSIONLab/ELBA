// Created by Saliya Ekanayake on 1/7/19.

#ifndef LBL_DAL_FASTADATA_HPP
#define LBL_DAL_FASTADATA_HPP

#include <cstdlib>
#include <memory>
#include <vector>
#include "Types.hpp"

/*!
 * Utility to store and retrive the data read from FASTA files.
 */
class FastaData {
public:
  /*!
   * Creates a new FastaData instance. Note, data could contain more
   * characters than required for this process. Therefore, l_start
   * and l_end identify the start and end (inclusive) of the data segment
   * relevant to this process. Note, l_end is expected to point to the last
   * character of the last sequence relevant to this process, so it doesn't
   * include a new line character at the end.
   *
   * While creating the FastData instance, any sequence with length less than
   * the k-mer size (k) is removed. Also, multi-line FASTA content is
   * modified in-place to be single lined. Moreover, any * characters that
   * can exist in the middle of a sequence is removed.
   *
   * The above removals will modify the l_end pointer as well the total number
   * of sequences used in the computation.
   * @param data Shared pointer to the raw character stream. Note, this may
   * contain more characters than necessary by the owning process.
   * @param l_start Starting offset of the data for this process.
   * @param l_end End offset (inclusive) of the data for this process.
   * @param k The k-mer size
   */
  FastaData(std::shared_ptr<char> data, uint64_t l_start, uint64_t l_end,
            ushort k);

  /*!
   * Destructor for FastaData
   */
  ~FastaData();

  /*!
   * Returns a shared pointer to the raw FASTA data stream, and sets the
   * requested sequence's length, starting offset, and end offset (inclusive)
   * in the output parameters.
   * @param idx Requested sequence index.
   * @param is_global Boolean flag to indicate if the requested index is global
   * (the original sequence index in the file) or local (sequences local to
   * this process only).
   * @param [out] len The length of the sequence.
   * @param [out] start_offset The start offset of the sequence
   * @param [out] end_offset_inclusive The end offset (inclusive) of
   * the sequence
   * @return A shared pointer to the raw character stream of the FASTA content.
   */
  std::shared_ptr<char> get_sequence(uint64_t idx, bool is_global, ushort &len,
                                     uint64_t &start_offset,
                                     uint64_t &end_offset_inclusive);

  /*!
   *
   * @return The number of sequences local to this process.
   */
  uint64_t count();

  /*!
   *
   * @return The global sequence offset
   */
  uint64_t offset();

  /*!
   *
   * @return The global sequence count after removing any that is
   * less than the k-mer length
   */
  uint64_t global_count();

  /*!
   * A helper method to print information of the FastaData instance.
   */
  void print();

private:
  /*! Shared pointer to the raw character stream of FASTA content  */
  std::shared_ptr<char> data;
  /*! Number of sequences local to this process available in data stream */
  uint64_t l_seq_count;
  /*!
   * Global sequence count, which may be different from the total input sequence
   * count because some might get removed if their lengths are less than the
   * k-mer length.
   */
  uint64_t g_seq_count;
  /*! The global sequence offset, i.e. the sequence index of the
   * original FASTA file corresponding to the first sequence stored in this
   * instance.
   */
  uint64_t g_seq_offset;
  /*! Starting and ending offsets of the sequence data local to this instance
   * in the data stream
   */
  uint64_t l_start, l_end;
  /*! Offsets in the raw data stream to sequence ID starts */
  uvec_64* id_starts = nullptr;
  /*! Offsets in the raw data stream to the sequence (actual bases) starts */
  uvec_64* seq_starts = nullptr;


};


#endif //LBL_DAL_FASTADATA_HPP
