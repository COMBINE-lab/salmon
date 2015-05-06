#ifndef LINEARSYSTEMBUILDER_HPP
#define LINEARSYSTEMBUILDER_HPP

#include <cstdint>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <limits>
#include <iomanip>
#include <sstream>

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>

/** This class builds the linear system that represents the
 * transcript abundance problem.  Given the transcript database (T)
 * and the read database (R), let K|R be the set of kmers occurring
 * in R.  We then construct two quantities, A, a matrix and b, a vector.
 * Let N = \card{K|R} and M = \card{T}, then A is an NxM matrix and b is
 * an Nx1 vector.  Eventually, we will solve this system directly, but, for now,
 * we'll write it to file so we can call an external solver on it.
 */
template <typename HashT>
class LinearSystemBuilder {

  private:
  std::string _transcriptFile;
  HashT _readHash;
  std::unordered_map<uint64_t, uint64_t> _indexMap;
  size_t _numTranscripts;

  bool buildIndexMap() {
    uint_t merLen = _readHash.get_mer_len();

    const char* fnames[] = { _transcriptFile.c_str() };
    jellyfish::parse_read parser( fnames, fnames+1, 100);
    jellyfish::parse_read::thread stream = parser.new_thread();

    // Read pointer
    jellyfish::read_parser::read_t *read;
    // Read the input data
    while ( (read = stream.next_read()) ) {

      // Ignore header for now
      //std::string header( read->header, read->hlen );
      std::string seq(read->seq_s, std::distance(read->seq_s, read->seq_e)-1 );

      // Remove newlines from the read
      auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
      auto readLen = std::distance( seq.begin(), newEnd );

      size_t offset = 0;
      size_t mapped = 0;
      uint64_t indexSize = _indexMap.size();

      while ( offset < readLen - merLen )  {
        auto mer = seq.substr( offset, merLen );
        auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
        // add it if we don't find it
        if ( _indexMap.find(binMer) == _indexMap.end() ) {
          _indexMap[binMer] = indexSize;
          // We've grown the index
          ++indexSize;
        }
        ++offset;
      } // end kmer loop

      std::cerr << "processed " << _numTranscripts << " transcripts\n";

      ++_numTranscripts;
    } // end transcripts loop

    return true;
  }

  public:
  LinearSystemBuilder( std::string& transcriptFile, HashT& readHash) :
    _transcriptFile( transcriptFile ), _readHash( readHash ), _numTranscripts( 0 ) {}

  void writeToFiles( const std::string& matFileName, const std::string& bFileName) {
    // Map each occurring kmer to an index
    buildIndexMap();

    auto numRows = _indexMap.size();
    auto numCols = _numTranscripts;

    {
      // Put the right-hand side in a vector since we'll have
      // to write it out in order
      std::vector<float> b(numRows,0.0f);
      for ( auto& rowInd : _indexMap ){ b[rowInd.second] = static_cast<float>(_readHash[rowInd.first]); }

      // write it out
      std::ofstream bFile(bFileName);
      for ( auto be : b ) { bFile << be << "\n"; }
      bFile.close();
    }

    uint_t merLen = _readHash.get_mer_len();
    const char* fnames[] = { _transcriptFile.c_str() };
    jellyfish::parse_read parser( fnames, fnames+1, 100);
    jellyfish::parse_read::thread stream = parser.new_thread();

    std::ofstream matFile(matFileName);

    std::stringstream ss(std::stringstream::out);
    ss << "SPARSE\n";
    ss << numRows << " " << numCols <<  "\n";

    matFile << ss.str();

    size_t nnzFileOffset = ss.str().length();

    ss.clear();
    ss << std::numeric_limits<size_t>::max();
    auto nnzLen = ss.str().length();
    std::string placeHolder(nnzLen, ' ');
    matFile << "\n";
    //matFile << placeHolder << "\n";

    size_t nnz{0};
    size_t col{0};
    // Read pointer
    jellyfish::read_parser::read_t *read;
    // Read the input data
    while ( (read = stream.next_read()) ) {

      // Ignore header for now
      std::string header( read->header, read->hlen );
      std::string seq(read->seq_s, std::distance(read->seq_s, read->seq_e)-1 );

      // Remove newlines from the read
      auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
      auto readLen = std::distance( seq.begin(), newEnd );

      size_t offset = 0;
      size_t mapped = 0;
      uint64_t indexSize = _indexMap.size();

      std::unordered_map<size_t, size_t> rowEntries;
      while ( offset < readLen - merLen )  {
        auto mer = seq.substr( offset, merLen );
        auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), merLen );
        auto row = _indexMap[binMer];
        rowEntries[row] += 1;
        ++offset;
      } // end kmer loop

      // write out entries
      for ( auto& rowVal : rowEntries ) {
        matFile << col+1 << "\t" << rowVal.first+1 << "\t" << rowVal.second << "\n";
        ++nnz;
      }
      std::cerr << "processed " << col << " columns\n";
      ++col;

    } // end transcripts loop

    matFile.seekp(nnzFileOffset);
    matFile << nnz << "\n";
    matFile.close();
  }


};

#endif //LINEARSYSTEMBUILDER_HPP
