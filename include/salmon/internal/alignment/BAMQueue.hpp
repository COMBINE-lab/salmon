#ifndef __BAMQUEUE_HPP__
#define __BAMQUEUE_HPP__

//#include <boost/lockfree/spsc_queue.hpp>
//#include <boost/lockfree/queue.hpp>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

#include "salmon/internal/alignment/AlignmentGroup.hpp"
#include "salmon/internal/model/LibraryFormat.hpp"
#include "salmon/internal/alignment/ReadPair.hpp"
#include "salmon/internal/util/SalmonMath.hpp"
#include "salmon/internal/alignment/UnpairedRead.hpp"
#include "salmon/vendor/concurrentqueue.h"
#include "salmon/vendor/readerwriterqueue.h"
#include <spdlog/spdlog.h>
#include <boost/filesystem.hpp>
#include <boost/timer/timer.hpp>
#include <oneapi/tbb/concurrent_queue.h>

#include "salmon/internal/io/AlignmentIO.hpp"

/**
 * Simple structure holding info about the alignment file.
 */
struct AlignmentFile {
  boost::filesystem::path fileName;
  std::string readMode;
  AlignmentFileHandle* fp;
  AlignmentHeader* header;
  uint32_t numParseThreads;
};

/**
 * A queue from which to draw BAM alignments.  The queue is thread-safe, and
 * can be written to and read from multiple threads.
 *
 * This class is templated on LibT --- the type of read library from which
 * the provided alignments are being generated.
 */
template <typename FragT> class BAMQueue {
public:
  BAMQueue(std::vector<boost::filesystem::path>& fnames, LibraryFormat& libFmt,
           uint32_t numParseThreads, uint32_t cacheSize);
  ~BAMQueue();
  void forceEndParsing();

  AlignmentHeader* header();
  AlignmentHeader* safeHeader();

  std::vector<AlignmentHeader*> headers();

  template <typename FilterT>
  void start(FilterT filt, bool onlyProcessAmbiguousAlignments = false);

  inline bool getAlignmentGroup(AlignmentGroup<FragT*>*& group);

  // Return the number of reads processed so far by the queue
  size_t numObservedAlignments();
  size_t numObservedFragments();
  size_t numMappedFragments();
  size_t numUniquelyMappedFragments();

  void reset();

  oneapi::tbb::concurrent_queue<FragT*>& getFragmentQueue();
  // moodycamel::ConcurrentQueue<FragT*>& getFragmentQueue();

  // oneapi::tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*>&
  // getAlignmentGroupQueue();
  moodycamel::ConcurrentQueue<AlignmentGroup<FragT*>*>&
  getAlignmentGroupQueue();

private:
  size_t popNum{0};
  /** Fill the queue with the appropriate type of alignment
   * depending on the template paramater T
   */
  template <typename FilterT> void fillQueue_(FilterT, bool);

  /** Overload of getFrag_ for paired-end reads */
  template <typename FilterT>
  inline bool getFrag_(ReadPair& rpair, FilterT filt);
  /** Overload of getFrag_ for single-end reads */
  template <typename FilterT>
  inline bool getFrag_(UnpairedRead& sread, FilterT filt);

  // ---- mate-interleaving read-ahead buffer (paired-end) ------------------
  // Some aligners (e.g. `bowtie2 -k/-a`) emit, within a name-collated fragment,
  // *all* read-1 records followed by *all* read-2 records rather than
  // interleaving each reported pair's two mates. The paired-end getFrag_ relies
  // on mates being adjacent in the stream, so without help it would pair
  // read1-with-read1. These helpers buffer one whole name group and reorder it
  // so each reported pair's mates (matched by RNEXT/PNEXT, disambiguated by HI)
  // are adjacent — a no-op when the input is already interleaved.

  /** Read the next raw record from the input, transparently advancing across
   *  input files; returns false only when all files are exhausted. */
  bool rawNextRecord_(bam_seq_t** dest);
  /** Buffer the next whole name-collated group and reorder it so mates are
   *  adjacent. Leaves `groupCount_ == 0` at end of input. */
  void fillGroupBuffer_();
  /** Serve the next record of the (mate-interleaved) group buffer, refilling
   *  as needed; returns false at end of input. */
  bool bufferedNextRecord_(bam_seq_t** dest);
  /** Discard buffered state (between parsing passes). */
  void resetGroupBuffer_();

  // owned, reusable record allocations for the current group
  std::vector<bam_seq_t*> recPool_;
  // serve order into recPool_ (a permutation of [0, groupCount_) with mates adjacent)
  std::vector<size_t> groupOrder_;
  size_t groupCount_{0}; // valid records in recPool_ for the current group
  size_t groupPos_{0};   // next index into groupOrder_ to serve
  bam_seq_t* carryRec_{nullptr}; // one-record look-ahead (first record of next group)
  bool haveCarry_{false};
  std::vector<char> groupUsed_; // scratch: which records have been paired

public:
  bool verbose = false;

private:
  std::vector<AlignmentFile> files_;
  std::string fname_;
  LibraryFormat libFmt_;

  std::vector<AlignmentFile>::iterator currFile_;
  AlignmentFileHandle* fp_ = nullptr;
  AlignmentHeader* hdr_ = nullptr;

  // htsFile* fp_ = nullptr;
  size_t totalAlignments_;
  size_t numUnaligned_;
  size_t numMappedReads_;
  size_t numUniquelyMappedReads_;
  oneapi::tbb::concurrent_queue<FragT*> fragmentQueue_;
  // moodycamel::ConcurrentQueue<FragT*> fragmentQueue_;

  // oneapi::tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*> alnGroupPool_;
  moodycamel::ConcurrentQueue<AlignmentGroup<FragT*>*> alnGroupPool_;

  // oneapi::tbb::concurrent_bounded_queue<AlignmentGroup<FragT*>*> alnGroupQueue_;
  moodycamel::ReaderWriterQueue<AlignmentGroup<FragT*>*> alnGroupQueue_;

  /*
  boost::lockfree::spsc_queue<AlignmentGroup<FragT*>*,
                              boost::lockfree::capacity<65535>> alnGroupQueue_;
                              */
  volatile bool doneParsing_;
  volatile bool exhaustedAlnGroupPool_;
  std::unique_ptr<std::thread> parsingThread_;
  std::shared_ptr<spdlog::logger> logger_;

  size_t batchNum_;
  std::string readMode_;
  /*
  #if not defined(__APPLE__)
     std::mutex agMutex_;
      std::condition_variable workAvailable_;
  #endif
  */
};

#include "salmon/internal/alignment/BAMQueue.tpp"
#endif //__BAMQUEUE_HPP__
