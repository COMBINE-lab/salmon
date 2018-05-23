#ifndef __FASTX_PARSER__
#define __FASTX_PARSER__

#include "fcntl.h"
#include "unistd.h"
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <thread>
#include <vector>

extern "C" {
#include "kseq.h"
}

#include "concurrentqueue.h"

#ifndef __FASTX_PARSER_PRECXX14_MAKE_UNIQUE__
#define __FASTX_PARSER_PRECXX14_MAKE_UNIQUE__

#if __cplusplus >= 201402L
#include <memory>
using std::make_unique;
#else

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

template <class T> struct _Unique_if {
  using _Single_object = std::unique_ptr<T>;
};

template <class T> struct _Unique_if<T[]> {
  using _Unknown_bound = std::unique_ptr<T[]>;
};

template <class T, size_t N> struct _Unique_if<T[N]> {
  using _Known_bound = void;
};

template <class T, class... Args>
typename _Unique_if<T>::_Single_object make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template <class T>
typename _Unique_if<T>::_Unknown_bound make_unique(size_t n) {
  using U = typename std::remove_extent<T>::type;
  return std::unique_ptr<T>(new U[n]());
}

template <class T, class... Args>
typename _Unique_if<T>::_Known_bound make_unique(Args&&...) = delete;

#endif // C++11
#endif //__FASTX_PARSER_PRECXX14_MAKE_UNIQUE__

namespace fastx_parser {
struct ReadSeq {
    std::string seq;
    std::string name;
    ~ReadSeq() {}
};

struct ReadQual {
  std::string seq;
  std::string name;
  std::string qual;
  ~ReadQual() {}
};

struct ReadPair {
  ReadSeq first;
  ReadSeq second;
};

struct ReadQualPair {
  ReadQual first;
  ReadQual second;
};

template <typename T> class ReadChunk {
public:
  ReadChunk(size_t want) : group_(want), want_(want), have_(want) {}
  inline void have(size_t num) { have_ = num; }
  inline size_t size() { return have_; }
  inline size_t want() const { return want_; }
  T& operator[](size_t i) { return group_[i]; }
  typename std::vector<T>::iterator begin() { return group_.begin(); }
  typename std::vector<T>::iterator end() { return group_.begin() + have_; }

private:
  std::vector<T> group_;
  size_t want_;
  size_t have_;
};

template <typename T> class ReadGroup {
public:
  ReadGroup(moodycamel::ProducerToken&& pt, moodycamel::ConsumerToken&& ct)
      : pt_(std::move(pt)), ct_(std::move(ct)) {}
  moodycamel::ConsumerToken& consumerToken() { return ct_; }
  moodycamel::ProducerToken& producerToken() { return pt_; }
  // get a reference to the chunk this ReadGroup owns
  std::unique_ptr<ReadChunk<T>>& chunkPtr() { return chunk_; }
  // get a *moveable* reference to the chunk this ReadGroup owns
  std::unique_ptr<ReadChunk<T>>&& takeChunkPtr() { return std::move(chunk_); }
  inline void have(size_t num) { chunk_->have(num); }
  inline size_t size() { return chunk_->size(); }
  inline size_t want() const { return chunk_->want(); }
  T& operator[](size_t i) { return (*chunk_)[i]; }
  typename std::vector<T>::iterator begin() { return chunk_->begin(); }
  typename std::vector<T>::iterator end() {
    return chunk_->begin() + chunk_->size();
  }
  void setChunkEmpty() { chunk_.release(); }
  bool empty() const { return chunk_.get() == nullptr; }

private:
  std::unique_ptr<ReadChunk<T>> chunk_{nullptr};
  moodycamel::ProducerToken pt_;
  moodycamel::ConsumerToken ct_;
};

template <typename T> class FastxParser {
public:
  FastxParser(std::vector<std::string> files, uint32_t numConsumers,
              uint32_t numParsers = 1, uint32_t chunkSize = 1000);

  FastxParser(std::vector<std::string> files, std::vector<std::string> files2,
              uint32_t numConsumers, uint32_t numParsers = 1,
              uint32_t chunkSize = 1000);
  ~FastxParser();
  bool start();
  bool stop();
  ReadGroup<T> getReadGroup();
  bool refill(ReadGroup<T>& rg);
  void finishedWithGroup(ReadGroup<T>& s);

private:
  moodycamel::ProducerToken getProducerToken_();
  moodycamel::ConsumerToken getConsumerToken_();

  std::vector<std::string> inputStreams_;
  std::vector<std::string> inputStreams2_;
  uint32_t numParsers_;
  std::atomic<uint32_t> numParsing_;

  // NOTE: Would like to use std::future<int> here instead, but that
  // solution doesn't seem to work.  It's unclear exactly why
  // see (https://twitter.com/nomad421/status/917748383321817088)
  std::vector<std::unique_ptr<std::thread>> parsingThreads_;

  // holds the results of the parsing threads, which is simply equal to
  // the return value of kseq_read() for the last call to that function.
  // A value < -1 signifies some sort of error.
  std::vector<int> threadResults_;

  size_t blockSize_;
  moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>> readQueue_,
      seqContainerQueue_;

  // holds the indices of files (file-pairs) to be processed
  moodycamel::ConcurrentQueue<uint32_t> workQueue_;

  std::vector<std::unique_ptr<moodycamel::ProducerToken>> produceReads_;
  std::vector<std::unique_ptr<moodycamel::ConsumerToken>> consumeContainers_;
  bool isActive_{false};
};
}
#endif // __FASTX_PARSER__
