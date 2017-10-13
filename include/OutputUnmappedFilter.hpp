#ifndef __OUTPUT_UNMAPPED_FILTER_HPP__
#define __OUTPUT_UNMAPPED_FILTER_HPP__

#include "ReadPair.hpp"
#include "UnpairedRead.hpp"
#include <iostream>
#include <tbb/concurrent_queue.h>

template <typename FragT> class OutputUnmappedFilter {
public:
  OutputUnmappedFilter(tbb::concurrent_bounded_queue<FragT*>* outQueue)
      : outQueue_(outQueue), qlen_(0) {
    memset(&qname_[0], 0, 255);
  }

  inline void processFrag(FragT* f) {
    bool distinctRead{true};

    // Check if the new read name is distinct from
    // the previous one.
    auto newLen = f->getNameLength();
    if (qlen_ == newLen) {
      distinctRead = (memcmp(&qname_[0], f->getName(), qlen_) != 0);
    }

    // If this doesn't have the same name as the last read we saw,
    // than pass it to the output queue and update the name of the
    // most recently observed read.
    if (distinctRead) {
      auto* alnCpy = f->clone();
      outQueue_->push(alnCpy);
      qlen_ = newLen;
      memcpy(&qname_[0], f->getName(), qlen_);
    }
    // If it does have the same name as the last read, no need
    // to update.
  }

private:
  tbb::concurrent_bounded_queue<FragT*>* outQueue_ = nullptr;
  char qname_[255];
  uint32_t qlen_;
};

#endif // __OUTPUT_UNMAPPED_FILTER_HPP__
