/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef __JELLYFISH_COOPERATIVE_POOL_HPP__
#define __JELLYFISH_COOPERATIVE_POOL_HPP__

#include <assert.h>
#include <unistd.h>

#include <jellyfish/circular_buffer.hpp>
#include <jellyfish/compare_and_swap.hpp>
#include <jellyfish/locks_pthread.hpp>

/// Cooperative pool. Provide a link between a producer and many
/// consumers. It is cooperative in the sense that there is no
/// dedicated thread to the producer. When the number of elements in
/// the queue from the producer to the consumer is less than half,
/// then the thread requesting an element attempts to become the
/// producer. It stays a producer until the producer to consumer queue
/// is full.
///
/// This class must be subclassed using CRTP. `T` is the type of the
/// element passed around in the queues. The derived class must
/// implement the method `bool produce(T& e)`. It is called when a
/// thread has become a producer. It must set in `e` the new element,
/// unless there is nothing more to produce. It returns `true` if
/// there is nothing more to produce (and `e` is not used), `false`
/// otherwise.
///
/// The following example will produce the integers `[0, 1000000]`:
///
/// ~~~{.cc}
/// class sequence : public cooperative_bool<sequence, int> {
///   int cur_;
/// public:
///   sequence() : cur_(0) { }
///   bool produce(int& e) {
///     if(cur_ <= 1000000) {
///       e = cur_++;
///       return false;
///     }
///     return true;
///   }
/// };
/// ~~~
///
/// To access the elements (or the jobs) of the sequence, instantiate
/// a `sequence::job` object and check that it is not empty. If empty,
/// the sequence is over.
///
/// ~~~{.cc}
/// sequence seq; // Sequence, instantiated in main thread
/// // In each consumer thread:
/// while(true) {
///   sequence::job j(seq);
///   if(j.is_empty())
///     break;
///   // Do computation using *j and j->
/// }
/// ~~~

namespace jellyfish {
template<typename D, typename T>
class cooperative_pool {
public:
  typedef jflib::circular_buffer<uint32_t> cbT;
  typedef T                                element_type;

private:
  uint32_t      size_;
  element_type* elts_;
  cbT           cons_prod_;     // FIFO from Consumers to Producers
  cbT           prod_cons_;     // FIFO from Producers to Consumers
  int           has_producer_;  // Tell whether a thread is acting as a producer

  // RAII token.
  class take_token {
    int* const token_;
    const bool has_token_;
  public:
    take_token(int* token) : token_(token), has_token_(jflib::cas(token_, 0, 1)) { }
    ~take_token() {
      if(has_token_)
        //        cas(token_, 1, 0); // Guaranteed to succeed. Memory barrier
        jflib::a_store(token_, 0);
    }
    bool has_token() const { return has_token_; }
  };

  explicit cooperative_pool(const cooperative_pool& rhs) : size_(0), elts_(0), cons_prod_(0), prod_cons_(0), has_producer_(0) { }
public:
  cooperative_pool(uint32_t size) :
    size_(size),
    elts_(new element_type[size_]),
    cons_prod_(size_ + 100),
    prod_cons_(size_ + 100),
    has_producer_(0)
  {
    // Every element is empty and ready to be filled by the producer
    for(size_t i = 0; i < size_; ++i)
      cons_prod_.enqueue_no_check(i);
  }

  ~cooperative_pool() { delete [] elts_; }

  uint32_t size() const { return size_; }

  element_type* element_begin() { return elts_; }
  element_type* element_end() { return elts_ + size_; }

  // Contains a filled element or is empty. In which case the producer
  // is done and we should stop processing.
  class job {
    cooperative_pool& cp_;
    uint32_t          i_;       // Index of element
  public:
    job(cooperative_pool& cp) : cp_(cp), i_(cp_.get_element()) { }
    ~job() { release(); }

    void release() {
      if(!is_empty()) {
        cp_.cons_prod_.enqueue_no_check(i_);
      }
    }
    bool is_empty() const { return i_ == cbT::guard; }
    void next() {
      release();
      i_ = cp_.get_element();
    }

    element_type& operator*() { return cp_.elts_[i_]; }
    element_type* operator->() { return &cp_.elts_[i_]; }

  private:
    // Disable copy of job
    job(const job& rhs) { }
    job& operator=(const job& rhs) { }
  };
  friend class job;

  /// STL compliant iterator
  class iterator : public std::iterator<std::input_iterator_tag, element_type> {
    job* j_;
  public:
    iterator() : j_(0) { }
    iterator(cooperative_pool& cp) : j_(new job(cp)) { }
    iterator(const iterator& rhs) : j_(rhs.j_) { }

    bool operator==(const iterator& rhs) const { return j_ == rhs.j_; }
    bool operator!=(const iterator& rhs) const { return j_ != rhs.j_; }
    element_type& operator*() { return j_->operator*(); }
    element_type* operator->() { return j_->operator->(); }

    iterator& operator++() {
      j_->next();
      if(j_->is_empty()) {
        delete j_;
        j_ = 0;
      }
      return *this;
    }

    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  };
  iterator begin() { return iterator(*this); }
  const iterator begin() const { return iterator(*this); }
  const iterator end() const { return iterator(); }

private:
  enum PRODUCER_STATUS { PRODUCER_PRODUCED, PRODUCER_DONE, PRODUCER_EXISTS };
  uint32_t get_element() {
    int iteration = 0;

    while(true) {
      // If less than half full -> try to fill up producer to consumer
      // queue. Disregard return value: in any every case will
      // attempt to get an element for ourselves
      if(prod_cons_.fill() < prod_cons_.size() / 2)
        become_producer();

      uint32_t i = prod_cons_.dequeue();
      if(i != cbT::guard)
        return i;

      // Try to become producer
      switch(become_producer()) {
      case PRODUCER_PRODUCED:
        iteration = 0; // Produced. Attempt anew to get an element
        break;
      case PRODUCER_DONE:
        return prod_cons_.dequeue();
      case PRODUCER_EXISTS:
        delay(iteration++); // Already a producer. Wait a bit it adds things to queue
        break;
      }
    }
  }

  PRODUCER_STATUS become_producer() {
    if(prod_cons_.is_closed())
      return PRODUCER_DONE;

    // Mark that we have a produce (myself). If not, return. Token
    // will be release automatically at end of method.
    take_token producer_token(&has_producer_);
    if(!producer_token.has_token())
      return PRODUCER_EXISTS;

    uint32_t i = cbT::guard;
    try {
      while(true) { // Only way out is if produce method is done (returns true or throw an exception)
        i = cons_prod_.dequeue();
        if(i == cbT::guard)
          return PRODUCER_PRODUCED;

        if(static_cast<D*>(this)->produce(elts_[i])) // produce returns true if done
          break;

        prod_cons_.enqueue_no_check(i);
      }
    } catch(...) { }       // Threw an exception -> same as being done

    // Producing is done
    cons_prod_.enqueue_no_check(i);
    prod_cons_.close();

    return PRODUCER_DONE;
  }

  // First 16 operations -> no delay. Then exponential back-off up to a second.
  void delay(int iteration) {
    if(iteration < 16)
      return;
    int shift = 10 - std::min(iteration - 16, 10);
    usleep((1000000 - 1) >> shift);
  }
};

} // namespace jellyfish {
#endif /* __JELLYFISH_COOPERATIVE_POOL_HPP__ */
