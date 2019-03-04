#include "FastxParser.hpp"
#include "FastxParserThreadUtils.hpp"

#include "fcntl.h"
#include "unistd.h"
#include <sstream>
#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <poll.h>
#include <thread>
#include <vector>
#include <zlib.h>

// STEP 1: declare the type of file handler and the read() function
KSEQ_INIT(gzFile, gzread)

namespace fastx_parser {
template <typename T>
FastxParser<T>::FastxParser(std::vector<std::string> files,
                            uint32_t numConsumers, uint32_t numParsers,
                            uint32_t chunkSize)
    : FastxParser(files, {}, numConsumers, numParsers, chunkSize) {}

template <typename T>
FastxParser<T>::FastxParser(std::vector<std::string> files,
                            std::vector<std::string> files2,
                            uint32_t numConsumers, uint32_t numParsers,
                            uint32_t chunkSize)
    : inputStreams_(files), inputStreams2_(files2), numParsing_(0),
      blockSize_(chunkSize) {

  if (numParsers > files.size()) {
    std::cerr << "Can't make user of more parsing threads than file (pairs); "
                 "setting # of parsing threads to "
              << files.size();
    numParsers = files.size();
  }
  numParsers_ = numParsers;

  // nobody is parsing yet
  numParsing_ = 0;

  readQueue_ = moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>(
      4 * numConsumers, numParsers, 0);

  seqContainerQueue_ =
      moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>(
          4 * numConsumers, 1 + numConsumers, 0);

  workQueue_ = moodycamel::ConcurrentQueue<uint32_t>(numParsers_);

  // push all file ids on the queue
  for (size_t i = 0; i < files.size(); ++i) {
    workQueue_.enqueue(i);
  }

  // every parsing thread gets a consumer token for the seqContainerQueue
  // and a producer token for the readQueue.
  for (size_t i = 0; i < numParsers_; ++i) {
    consumeContainers_.emplace_back(
        new moodycamel::ConsumerToken(seqContainerQueue_));
    produceReads_.emplace_back(new moodycamel::ProducerToken(readQueue_));
  }

  // enqueue the appropriate number of read chunks so that we can start
  // filling them once the parser has been started.
  moodycamel::ProducerToken produceContainer(seqContainerQueue_);
  for (size_t i = 0; i < 4 * numConsumers; ++i) {
    auto chunk = make_unique<ReadChunk<T>>(blockSize_);
    seqContainerQueue_.enqueue(produceContainer, std::move(chunk));
  }
}

template <typename T> ReadGroup<T> FastxParser<T>::getReadGroup() {
  return ReadGroup<T>(getProducerToken_(), getConsumerToken_());
}

template <typename T>
moodycamel::ProducerToken FastxParser<T>::getProducerToken_() {
  return moodycamel::ProducerToken(seqContainerQueue_);
}

template <typename T>
moodycamel::ConsumerToken FastxParser<T>::getConsumerToken_() {
  return moodycamel::ConsumerToken(readQueue_);
}

template <typename T> FastxParser<T>::~FastxParser() {
  if (isActive_ or numParsing_ > 0) {
    // Think about if this is too noisy --- but the user really shouldn't do this.
    std::cerr << "\n\nEncountered FastxParser destructor while parser was still marked active (or while parsing threads were still active). "
              << "Be sure to call stop() before letting FastxParser leave scope!\n";
    try {
      stop();
    } catch (const std::exception& e) {
      // Should exiting here be a user-definable behavior?
      // What is the right mechanism for that.
      std::cerr << "\n\nParser encountered exception : " << e.what() << "\n";
      std::exit(-1);
    }
  }
  // Otherwise, we are good to go (i.e., destruct)
}

  template <typename T>
  bool FastxParser<T>::stop() {
    bool ret{false};
    if (isActive_) {
      for (auto& t : parsingThreads_) {
        t->join();
      }
      isActive_ = false;
      for (auto& res : threadResults_) {
        if (res == -3) {
          throw std::range_error("Error reading from the FASTA/Q stream. Make sure the file is valid.");
        } else if (res < -1) {
          std::stringstream ss;
          ss << "Error reading from the FASTA/Q stream. Minimum return code for left and right read was ("
             << res << "). Make sure the file is valid.";
          throw std::range_error(ss.str());
        }
      }
      ret = true;
    } else {
      // Is this being too loud?  Again, if this triggers, the user has violated the API.
      std::cerr << "stop() was called on a FastxParser that was not marked active. Did you remember "
                << "to call start() on this parser?\n";
    }
    return ret;
  }

inline void copyRecord(kseq_t* seq, ReadSeq* s) {
  // Copy over the sequence and read name
  s->seq.assign(seq->seq.s, seq->seq.l);
  s->name.assign(seq->name.s, seq->name.l);
}

inline void copyRecord(kseq_t* seq, ReadQual* s) {
    // Copy over the sequence and read name 
    // and quality
    s->seq.assign(seq->seq.s, seq->seq.l);
    s->name.assign(seq->name.s, seq->name.l);
    s->qual.assign(seq->qual.s, seq->qual.l);
}


template <typename T>
int parseReads(
    std::vector<std::string>& inputStreams, std::atomic<uint32_t>& numParsing,
    moodycamel::ConsumerToken* cCont, moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<uint32_t>& workQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>&
        seqContainerQueue_,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue_) {

  using fastx_parser::thread_utils::MIN_BACKOFF_ITERS;
  auto curMaxDelay = MIN_BACKOFF_ITERS;
  kseq_t* seq;
  T* s;
  uint32_t fn{0};
  while (workQueue.try_dequeue(fn)) {
    auto file = inputStreams[fn];
    std::unique_ptr<ReadChunk<T>> local;
    while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
      fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      // Think of a way to do this that wouldn't be loud (or would allow a user-definable logging mechanism)
      // std::cerr << "couldn't dequeue read chunk\n";
    }
    size_t numObtained{local->size()};
    // open the file and init the parser
    auto fp = gzopen(file.c_str(), "r");

    // The number of reads we have in the local vector
    size_t numWaiting{0};

    seq = kseq_init(fp);
    int ksv = kseq_read(seq);

    while (ksv >= 0) {
      s = &((*local)[numWaiting++]);

      copyRecord(seq, s);

      // If we've filled the local vector, then dump to the concurrent queue
      if (numWaiting == numObtained) {
        curMaxDelay = MIN_BACKOFF_ITERS;
        while (!readQueue_.try_enqueue(std::move(local))) {
          fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
        }
        numWaiting = 0;
        numObtained = 0;
        // And get more empty reads
        curMaxDelay = MIN_BACKOFF_ITERS;
        while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
          fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
        }
        numObtained = local->size();
      }
      ksv = kseq_read(seq);
    }

    if (ksv == -3) {
      --numParsing;
      return -3;
    } else if (ksv < -1) {
      --numParsing;
      return ksv;
    }

    // If we hit the end of the file and have any reads in our local buffer
    // then dump them here.
    if (numWaiting > 0) {
      local->have(numWaiting);
      curMaxDelay = MIN_BACKOFF_ITERS;
      while (!readQueue_.try_enqueue(*pRead, std::move(local))) {
        fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      }
      numWaiting = 0;
    } else if (numObtained > 0){
      curMaxDelay = MIN_BACKOFF_ITERS;
      while (!seqContainerQueue_.try_enqueue(std::move(local))) {
        fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      }
    }
    // destroy the parser and close the file
    kseq_destroy(seq);
    gzclose(fp);
  }

  --numParsing;
  return 0;
}

template <typename T>
int parseReadPair(
    std::vector<std::string>& inputStreams,
    std::vector<std::string>& inputStreams2, std::atomic<uint32_t>& numParsing,
    moodycamel::ConsumerToken* cCont, moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<uint32_t>& workQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>&
        seqContainerQueue_,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue_) {

  using fastx_parser::thread_utils::MIN_BACKOFF_ITERS;
  size_t curMaxDelay = MIN_BACKOFF_ITERS;
  kseq_t* seq;
  kseq_t* seq2;
  T* s;

  uint32_t fn{0};
  while (workQueue.try_dequeue(fn)) {
    // for (size_t fn = 0; fn < inputStreams.size(); ++fn) {
    auto& file = inputStreams[fn];
    auto& file2 = inputStreams2[fn];

    std::unique_ptr<ReadChunk<T>> local;
    while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
      fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      // Think of a way to do this that wouldn't be loud (or would allow a user-definable logging mechanism)
      // std::cerr << "couldn't dequeue read chunk\n";
    }
    size_t numObtained{local->size()};
    // open the file and init the parser
    auto fp = gzopen(file.c_str(), "r");
    auto fp2 = gzopen(file2.c_str(), "r");

    // The number of reads we have in the local vector
    size_t numWaiting{0};

    seq = kseq_init(fp);
    seq2 = kseq_init(fp2);

    int ksv = kseq_read(seq);
    int ksv2 = kseq_read(seq2);
    while (ksv >= 0 and ksv2 >= 0) {

      s = &((*local)[numWaiting++]);
      copyRecord(seq, &s->first);
      copyRecord(seq2, &s->second);

      // If we've filled the local vector, then dump to the concurrent queue
      if (numWaiting == numObtained) {
        curMaxDelay = MIN_BACKOFF_ITERS;
        while (!readQueue_.try_enqueue(std::move(local))) {
          fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
        }
        numWaiting = 0;
        numObtained = 0;
        // And get more empty reads
        curMaxDelay = MIN_BACKOFF_ITERS;
        while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
          fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
        }
        numObtained = local->size();
      }
      ksv = kseq_read(seq);
      ksv2 = kseq_read(seq2);
    }

    if (ksv == -3 or ksv2 == -3) {
      --numParsing;
      return -3;
    } else if (ksv < -1 or ksv2 < -1) {
      --numParsing;
      return std::min(ksv, ksv2);
    }

    // If we hit the end of the file and have any reads in our local buffer
    // then dump them here.
    if (numWaiting > 0) {
      local->have(numWaiting);
      curMaxDelay = MIN_BACKOFF_ITERS;
      while (!readQueue_.try_enqueue(*pRead, std::move(local))) {
        fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      }
      numWaiting = 0;
    } else if (numObtained > 0){
      curMaxDelay = MIN_BACKOFF_ITERS;
      while (!seqContainerQueue_.try_enqueue(std::move(local))) {
        fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
      }
    }
    // destroy the parser and close the file
    kseq_destroy(seq);
    gzclose(fp);
    kseq_destroy(seq2);
    gzclose(fp2);
  }

  --numParsing;
  return 0;
}

template <> bool FastxParser<ReadSeq>::start() {
  if (numParsing_ == 0) {
    isActive_ = true;
    threadResults_.resize(numParsers_);
    std::fill(threadResults_.begin(), threadResults_.end(), 0);
    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
        this->threadResults_[i] = parseReads(this->inputStreams_, this->numParsing_,
                   this->consumeContainers_[i].get(),
                   this->produceReads_[i].get(), this->workQueue_,
                   this->seqContainerQueue_, this->readQueue_);
      }));
    }
    return true;
  } else {
    return false;
  }
}

template <> bool FastxParser<ReadPair>::start() {
  if (numParsing_ == 0) {
    isActive_ = true;
    // Some basic checking to ensure the read files look "sane".
    if (inputStreams_.size() != inputStreams2_.size()) {
      throw std::invalid_argument("There should be the same number "
                                  "of files for the left and right reads");
    }
    for (size_t i = 0; i < inputStreams_.size(); ++i) {
      auto& s1 = inputStreams_[i];
      auto& s2 = inputStreams2_[i];
      if (s1 == s2) {
        throw std::invalid_argument("You provided the same file " + s1 +
                                    " as both a left and right file");
      }
    }

    threadResults_.resize(numParsers_);
    std::fill(threadResults_.begin(), threadResults_.end(), 0);

    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
            this->threadResults_[i] = parseReadPair(this->inputStreams_, this->inputStreams2_,
                      this->numParsing_, this->consumeContainers_[i].get(),
                      this->produceReads_[i].get(), this->workQueue_,
                      this->seqContainerQueue_, this->readQueue_);
      }));
    }
    return true;
  } else {
    return false;
  }
}

template <> bool FastxParser<ReadQual>::start() {
    if (numParsing_ == 0) {
    isActive_ = true;
    threadResults_.resize(numParsers_);
    std::fill(threadResults_.begin(), threadResults_.end(), 0);
    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
        this->threadResults_[i] = parseReads(this->inputStreams_, this->numParsing_,
                   this->consumeContainers_[i].get(),
                   this->produceReads_[i].get(), this->workQueue_,
                   this->seqContainerQueue_, this->readQueue_);
      }));
    }
    return true;
  } else {
    return false;
  }
}

template <> bool FastxParser<ReadQualPair>::start() {
  if (numParsing_ == 0) {
    isActive_ = true;
    // Some basic checking to ensure the read files look "sane".
    if (inputStreams_.size() != inputStreams2_.size()) {
      throw std::invalid_argument("There should be the same number "
                                  "of files for the left and right reads");
    }
    for (size_t i = 0; i < inputStreams_.size(); ++i) {
      auto& s1 = inputStreams_[i];
      auto& s2 = inputStreams2_[i];
      if (s1 == s2) {
        throw std::invalid_argument("You provided the same file " + s1 +
                                    " as both a left and right file");
      }
    }

    threadResults_.resize(numParsers_);
    std::fill(threadResults_.begin(), threadResults_.end(), 0);

    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
            this->threadResults_[i] = parseReadPair(this->inputStreams_, this->inputStreams2_,
                      this->numParsing_, this->consumeContainers_[i].get(),
                      this->produceReads_[i].get(), this->workQueue_,
                      this->seqContainerQueue_, this->readQueue_);
      }));
    }
    return true;
  } else {
    return false;
  }
}


template <typename T> bool FastxParser<T>::refill(ReadGroup<T>& seqs) {
  finishedWithGroup(seqs);
  auto curMaxDelay = fastx_parser::thread_utils::MIN_BACKOFF_ITERS;
  while (numParsing_ > 0) {
    if (readQueue_.try_dequeue(seqs.consumerToken(), seqs.chunkPtr())) {
      return true;
    }
    fastx_parser::thread_utils::backoffOrYield(curMaxDelay);
  }
  return readQueue_.try_dequeue(seqs.consumerToken(), seqs.chunkPtr());
}

template <typename T> void FastxParser<T>::finishedWithGroup(ReadGroup<T>& s) {
  // If this read group is holding a valid chunk, then give it back
  if (!s.empty()) {
    seqContainerQueue_.enqueue(s.producerToken(), std::move(s.takeChunkPtr()));
    s.setChunkEmpty();
  }
}

template class FastxParser<ReadSeq>;
template class FastxParser<ReadPair>;
template class FastxParser<ReadQual>;
template class FastxParser<ReadQualPair>;
}
