#include "FastxParser.hpp"

#include "fcntl.h"
#include "unistd.h"
#include <atomic>
#include <cstdio>
#include <cstdlib>
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
  for (auto& t : parsingThreads_) {
    t->join();
  }
}

inline void copyRecord(kseq_t* seq, ReadSeq* s) {
  // Copy over the sequence and read name
  s->seq.assign(seq->seq.s, seq->seq.l);
  s->name.assign(seq->name.s, seq->name.l);
}

template <typename T>
void parseReads(
    std::vector<std::string>& inputStreams, std::atomic<uint32_t>& numParsing,
    moodycamel::ConsumerToken* cCont, moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<uint32_t>& workQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>&
        seqContainerQueue_,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue_) {
  kseq_t* seq;
  T* s;
  uint32_t fn{0};
  while (workQueue.try_dequeue(fn)) {
    auto file = inputStreams[fn];
    std::unique_ptr<ReadChunk<T>> local;
    while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
      std::cerr << "couldn't dequeue read chunk\n";
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
        while (!readQueue_.try_enqueue(std::move(local))) {
        }
        numWaiting = 0;
        numObtained = 0;
        // And get more empty reads
        while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
        }
        numObtained = local->size();
      }
      ksv = kseq_read(seq);
    }

    // If we hit the end of the file and have any reads in our local buffer
    // then dump them here.
    if (numWaiting > 0) {
      local->have(numWaiting);
      while (!readQueue_.try_enqueue(*pRead, std::move(local))) {
      }
      numWaiting = 0;
    }
    // destroy the parser and close the file
    kseq_destroy(seq);
    gzclose(fp);
  }

  --numParsing;
}

template <typename T>
void parseReadPair(
    std::vector<std::string>& inputStreams,
    std::vector<std::string>& inputStreams2, std::atomic<uint32_t>& numParsing,
    moodycamel::ConsumerToken* cCont, moodycamel::ProducerToken* pRead,
    moodycamel::ConcurrentQueue<uint32_t>& workQueue,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>&
        seqContainerQueue_,
    moodycamel::ConcurrentQueue<std::unique_ptr<ReadChunk<T>>>& readQueue_) {

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
      std::cerr << "couldn't dequeue read chunk\n";
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
        while (!readQueue_.try_enqueue(std::move(local))) {
        }
        numWaiting = 0;
        numObtained = 0;
        // And get more empty reads
        while (!seqContainerQueue_.try_dequeue(*cCont, local)) {
        }
        numObtained = local->size();
      }
      ksv = kseq_read(seq);
      ksv2 = kseq_read(seq2);
    }

    // If we hit the end of the file and have any reads in our local buffer
    // then dump them here.
    if (numWaiting > 0) {
      local->have(numWaiting);
      while (!readQueue_.try_enqueue(*pRead, std::move(local))) {
      }
      numWaiting = 0;
    }
    // destroy the parser and close the file
    kseq_destroy(seq);
    gzclose(fp);
    kseq_destroy(seq2);
    gzclose(fp2);
  }

  --numParsing;
}

template <> bool FastxParser<ReadSeq>::start() {
  if (numParsing_ == 0) {
    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
        parseReads(this->inputStreams_, this->numParsing_,
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
    for (size_t i = 0; i < numParsers_; ++i) {
      ++numParsing_;
      parsingThreads_.emplace_back(new std::thread([this, i]() {
        parseReadPair(this->inputStreams_, this->inputStreams2_,
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
  while (numParsing_ > 0) {
    if (readQueue_.try_dequeue(seqs.consumerToken(), seqs.chunkPtr())) {
      return true;
    }
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
}
