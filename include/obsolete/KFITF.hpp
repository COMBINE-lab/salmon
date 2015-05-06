#ifndef KFITF_HPP
#define KFITF_HPP

#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <mutex>
#include <unordered_set>

#include <jellyfish/sequence_parser.hpp>
#include <jellyfish/parse_read.hpp>
#include <jellyfish/mer_counting.hpp>
#include <jellyfish/misc.hpp>
#include <jellyfish/compacted_hash.hpp>


/** Kmer Frequence Inverse Transcript Frequency (KFITF) calculator
 *
 * We consider each k-mer as a "term" and each transcript as a "document".
 * The purpose of this class is to compute the TF-IDF (Term Frequency-Inverse Document Frequency)
 * of each kmer transcript pair.  Since the kmer frequency (KF) is on a per transcript basis,
 * it can be computed independently in each transcript and need not be stored in memory
 * after this transcript has been processed.  The inverse transcript freuqency (ITF), however
 * requires knowing the number of different transcripts in which a k-mer is located.  This can
 * be computed in a single pass over the data, but must be stored in memory.
 *
 */


class KFITFCalculator {
  private:
  typedef jellyfish::invertible_hash::array<uint64_t, atomic::gcc, allocators::mmap> hash_array;

  size_t _merLen;
  size_t _numTranscripts;
  hash_array _hash;

  public:
  KFITFCalculator( size_t merLen ) :
    _merLen( merLen ),
    _numTranscripts( 0 ),
    _hash( hash_array(1000000, // Size of hash. Will be rounded up to a power of 2
                     2*_merLen,      // Key length in bits (2 * mer length)
                     30,       // Value length in bits
                     126,     // Max # of reprobe
                     jellyfish::quadratic_reprobes // Reprobe policy
                     ) ) {

    //Test out the hash
    std::string mer{"AAAAAATTGTGTATGGC"};
    auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), _merLen );
    for ( size_t i = 0; i < 100; ++i ) {
      auto k = rand() % 100000 + 1;
      for ( size_t j = 0; j < k; ++j ) {
        _hash.add(binMer,1);
      }
      uint64_t val{0}; _hash.get_val(binMer, val);
      if( val != k ) {
        std::cerr << "val  = " << val << " but it should be  " << k  << "\n";
        std::abort();
      }
      _hash.map(binMer,0);
    }
  }

  template <typename StreamT>
  void computeITF( StreamT& stream ) {

    using std::string;

    std::mutex resLock;

    size_t numActors = 26;
    // Thread pool to count membership of multiple transcripts simultaneously
    boost::threadpool::pool tp(numActors);

    // Lock-free queue to hold kmer membership sets
    boost::lockfree::fifo< std::unordered_set<uint64_t>* > q;

    std::atomic<int> counter{0};

    auto readResults = [this, &tp, &counter] () -> bool {

      while ( tp.active() + tp.pending() > 1 ) {
        std::cerr << counter.load() << "\n";
        boost::this_thread::sleep_for(boost::chrono::milliseconds(250));
      }
      return true;
    };

    // schedule the consumer
    tp.schedule(readResults);


    // Read pointer
    jellyfish::read_parser::read_t *read;
    // Read the input data
    while ( (read = stream.next_read()) ) {
      std::string oseq(read->seq_s, std::distance(read->seq_s, read->seq_e)-1 );
      std::string header( read->header, read->hlen );

      auto task = [&resLock, &q, &counter, oseq, header, this]() -> void {
        std::string seq = oseq;
        auto newEnd  = std::remove( seq.begin(), seq.end(), '\n' );
        auto readLen = std::distance( seq.begin(), newEnd );

        std::unordered_map<uint64_t, size_t> occurringKmers;// = new std::unordered_set<uint64_t>;
        auto splitPos = header.find('|');
        auto geneName = header.substr(0, splitPos);
        size_t offset = 0;
        size_t mapped = 0;

        while ( offset < readLen - this->_merLen )  {
          auto mer = seq.substr( offset, this->_merLen );
          auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), this->_merLen );
          // add this kmer to the set for this transcript
          occurringKmers[binMer] += 1;//.insert( binMer );
          ++offset;
        }

        //resLock.lock();
        for ( const auto& kmer : occurringKmers ) {
          // Increment the count of this kmer
          this->_hash.add(kmer.first, kmer.second);
        }
        //resLock.unlock();

        ++counter;
        // Enqueue this transcript's kmer set
        //resLock.lock();
        //q.enqueue(occurringKmers);
        //resLock.unlock();
      }; // end task closure

      tp.schedule(task);

    } // end reading input

    // Wait for the work to finish
    tp.wait();
    _numTranscripts = counter.load();
  }

  float itf( const std::string& mer ) {
    auto binMer = jellyfish::parse_dna::mer_string_to_binary( mer.c_str(), _merLen );
    uint64_t val{0}; _hash.get_val( binMer, val );
    return 1.0 / (1.0f+val);
    if( _numTranscripts < val ) {
      std::cerr << "ERROR! Num Transcripts = " << _numTranscripts << ", but df(t,D) = " << val << "\n";
      std::abort();
    }
    return std::log2( static_cast<float>(_numTranscripts) / (1.0f + val));
  }

};

#endif //KFITF_HPP
