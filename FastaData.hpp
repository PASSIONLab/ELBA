// Created by Saliya Ekanayake on 1/7/19.

#ifndef LBL_DAL_FASTADATA_HPP
#define LBL_DAL_FASTADATA_HPP

#include <cstdlib>
#include <memory>
#include <vector>

/*!
 * Utility to store and retrive the data read from FASTA files.
 */
class FastaData {
public:
  /*!
   * Creates a new FastaData instance.
   * @param data Shared pointer to the raw character stream. Note, this may
   * contain more characters than necessary by the owning process.
   * @param l_start Starting offset of the data for this process.
   * @param l_end End offset (inclusive) of the data for this process.
   */
  FastaData(std::shared_ptr<char> data, int l_start, int l_end);

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
  std::shared_ptr<char> get_sequence(int idx, bool is_global, int &len,
                                     int &start_offset,
                                     int &end_offset_inclusive);

  /*!
   * A helper method to print information of the FastaData instance.
   */
  void print();

private:
  /*! Shared pointer to the raw character stream of FASTA content  */
  std::shared_ptr<char> data;
  /*! Number of sequences local to this process available in data stream */
  int l_seq_count;
  /*! The global sequence offset, i.e. the sequence index of the
   * original FASTA file corresponding to the first sequence stored in this
   * instance.
   */
  int g_seq_offset;
  /*! Starting and ending offsets of the sequence data local to this instance
   * in the data stream
   */
  int l_start, l_end;
  /*! Offsets in the raw data stream to sequence ID starts */
  std::vector<int>* id_starts = nullptr;
  /*! Offsets in the raw data stream to the sequence (actual bases) starts */
  std::vector<int>* seq_starts = nullptr;


};


#endif //LBL_DAL_FASTADATA_HPP
