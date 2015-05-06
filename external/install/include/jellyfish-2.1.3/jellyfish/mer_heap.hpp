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

#ifndef __JELLYFISH_HEAP_HPP__
#define __JELLYFISH_HEAP_HPP__

#include <stdint.h>
#include <algorithm>

namespace jellyfish { namespace mer_heap {
template<typename Key, typename Iterator>
struct heap_item {
  Key       key_;
  uint64_t  val_;
  uint64_t  pos_;
  Iterator* it_;

  heap_item() : it_(0) { }
  heap_item(Iterator& iter) : key_(iter.key()), val_(iter.val()), pos_(iter.pos()), it_(&iter) { }

  bool operator>(const heap_item& other) const {
    if(pos_ == other.pos_)
      return key_ > other.key_;
    return pos_ > other.pos_;
  }
};

// STL make_heap creates a max heap. We want a min heap, so
// we use the > operator
template<typename Key, typename Iterator>
struct heap_item_comp {
  bool operator()(const heap_item<Key, Iterator>* i1, const heap_item<Key, Iterator>* i2) {
    return *i1 > *i2;
  }
};

template<typename Key, typename Iterator>
class heap {
  heap_item<Key, Iterator>*     storage_; // Storage of the elements
  heap_item<Key, Iterator>**    elts_; // Pointers to storage. Create a heap of pointers
  size_t                        capacity_; // Actual capacity
  size_t                        h_; // Head pointer
  heap_item_comp<Key, Iterator> comp_;

public:
  typedef const heap_item<Key, Iterator> *const_item_t;

  heap() : storage_(0), elts_(0), capacity_(0), h_(0) { }
  explicit heap(size_t capacity)  { initialize(capacity); }
  ~heap() {
    delete[] storage_;
    delete[] elts_;
  }

  void initialize(size_t capacity) {
    capacity_ = capacity;
    h_        = 0;
    storage_   = new heap_item<Key, Iterator>[capacity_];
    elts_      = new heap_item<Key, Iterator>*[capacity_];
    for(size_t h1 = 0; h1 < capacity_; ++h1)
      elts_[h1] = &storage_[h1];
  }

  void fill(Iterator &it) {
    for(h_ = 0; h_ < capacity_; ++h_) {
      if(!it.next())
        break;
      storage_[h_] = it;
      elts_[h_]    = &storage_[h_];
    }
    std::make_heap(elts_, elts_ + h_, comp_);
  }
  // template<typename ForwardIterator>
  // void fill(ForwardIterator first, ForwardIterator last) {
  //   h_ = 0;
  //   while(h_ < capacity_ && first != last) {
  //     if(!first->next())
  //       break;
  //     storage_[h_].initialize(*first++);
  //     elts_[h_] = &storage_[h_];
  //     h_++;
  //   }
  //   std::make_heap(elts_, elts_ + h_, compare);
  // }

  bool is_empty() const { return h_ == 0; }
  bool is_not_empty() const { return h_ > 0; }
  size_t size() const { return h_; }
  size_t capacity() const { return capacity_; }

  // The following 3 should only be used after fill has been called
  const_item_t head() const { return elts_[0]; }
  void pop() { std::pop_heap(elts_, elts_ + h_--, comp_); }
  void push(Iterator &item) {
    *elts_[h_] = item;
    std::push_heap(elts_, elts_ + ++h_, comp_);
  }
};

} } // namespace jellyfish { namespace mer_heap {

#endif // __HEAP_HPP__
