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

#ifndef __LARGE_HASH_ITERATOR_HPP__
#define __LARGE_HASH_ITERATOR_HPP__

#include <iterator>
#include <utility>

/// Various iterators for the large hash array

namespace jellyfish { namespace large_hash {

/// Eager iterator. It computes the actual key and value when doing next.
template<typename array>
class eager_iterator_base {
public:
  typedef typename array::key_type    key_type;
  typedef typename array::mapped_type mapped_type;
  typedef typename array::key_status  key_status;

protected:
  const array* ary_;
  size_t       start_id_, id_, end_id_;
  key_type     key_;
  mapped_type  val_;

public:
  eager_iterator_base(const array* ary, size_t start, size_t end) :
    ary_(ary),
    start_id_(start > ary->size() ? ary->size() : start),
    id_(start),
    end_id_(end > ary->size() ? ary->size() : end)
  {}

  uint64_t start() const { return start_id_; }
  uint64_t end() const { return end_id_; }
  const key_type& key() const { return key_; }
  const mapped_type& val() const { return val_; }
  size_t id() const { return id_ - 1; }
  size_t pos() const { return key_.get_bits(0, ary_->lsize()); }

  bool next() {
    key_status success = array::EMPTY;
    while(success != array::FILLED && id_ < end_id_)
      success = ary_->get_key_val_at_id(id_++, key_, val_);
    if(success == array::FILLED)
      key_.set_bits(0, ary_->lsize(), ary_->inverse_matrix().times(key_));

    return success == array::FILLED;
  }
};

/// Lazy iterator. The actual key and value are actually computed only
/// when the key() and val() methods are called.
template<typename array>
class lazy_iterator_base {
public:
  typedef typename array::key_type    key_type;
  typedef typename array::mapped_type mapped_type;
  typedef typename array::key_status  key_status;
  typedef typename array::data_word   word;
  typedef typename array::offset_t    offset_t;

protected:
  const array*    ary_;
  size_t          start_id_, id_, end_id_;
  const word*     w_;
  const offset_t* o_;
  bool            reversed_key_;
  key_type        key_;

public:
  lazy_iterator_base(const array *ary, size_t start, size_t end) :
    ary_(ary),
    start_id_(ary ? (start > ary->size() ? ary->size() : start) : 0),
    id_(start),
    end_id_(ary ? (end > ary->size() ? ary->size() : end) : 0),
    w_(0), o_(0),
    reversed_key_(false)
  {}

  uint64_t start() const { return start_id_; }
  uint64_t end() const { return end_id_; }
  const key_type& key() {
    if(!reversed_key_) {
      key_.set_bits(0, ary_->lsize(), ary_->inverse_matrix().times(key_));
      reversed_key_ = true;
    }
    return key_;
  }
  mapped_type val() const {
    return ary_->get_val_at_id(id_ - 1, w_, o_, true, false);
  }
  size_t id() const { return id_ - 1; }
  size_t pos() const { return key_.get_bits(0, ary_->lsize()); }

  bool next() {
    reversed_key_      = false;
    key_status success = array::EMPTY;
    while(success != array::FILLED && id_ < end_id_)
      success = ary_->get_key_at_id(id_++, key_, &w_, &o_);

    return success == array::FILLED;
  }
};

/// Region iterator. Iterate over elements whose original position
/// (and not position after reprobing) falls inside the region
/// [start_id, end_id)
template<typename array>
class region_iterator_base {
  public:
  typedef typename array::key_type    key_type;
  typedef typename array::mapped_type mapped_type;
  typedef typename array::key_status  key_status;
  typedef typename array::data_word   word;
  typedef typename array::offset_t    offset_t;

protected:
  const array*    ary_;
  const uint64_t  mask_;
  const size_t    start_id_, end_id_, mid_;
  size_t          oid_, id_;
  const word*     w_;
  const offset_t* o_;
  bool            reversed_key_;
  key_type*       key_;
  bool            own_key;

public:
  region_iterator_base(const array *ary, size_t start, size_t end) :
    ary_(ary), mask_(ary ? ary->size() - 1 : 0),
    start_id_(ary ? std::min(start, ary->size()) : 0),
    end_id_(ary ? std::min(end, ary->size()) : 0),
    mid_(ary ? std::min(end_id_ - start_id_ + ary->max_reprobe_offset(), ary->size()) : 0),
    oid_(end_id_), id_(0), w_(0), o_(0),
    reversed_key_(false),
    key_(new key_type),
    own_key(true)
  {}

  region_iterator_base(const array *ary, size_t start, size_t end, key_type& key) :
    ary_(ary), mask_(ary ? ary->size() - 1 : 0),
    start_id_(ary ? std::min(start, ary->size()) : 0),
    end_id_(ary ? std::min(end, ary->size()) : 0),
    mid_(ary ? std::min(end_id_ - start_id_ + ary->max_reprobe_offset(), ary->size()) : 0),
    oid_(end_id_), id_(0), w_(0), o_(0),
    reversed_key_(false),
    key_(&key),
    own_key(false)
  { }

  ~region_iterator_base() {
    if(own_key)
      delete key_;
  }

  const key_type& key() {
    if(!reversed_key_) {
      key_->set_bits(0, ary_->lsize(), ary_->inverse_matrix().times(*key_));
      reversed_key_ = true;
    }
    return *key_;
  }
  mapped_type val() const {
    return ary_->get_val_at_id(id(), w_, o_, true, false);
  }
  uint64_t pos() const{
    return oid_;
  }

  size_t start() { return start_id_; }
  size_t end() { return end_id_; }

  /// Position where key is stored
  size_t id() const { return (start_id_ + id_ - 1) & mask_; }
  /// Original position (before reprobing).
  size_t oid() const { return oid_; }

  bool next() {
    reversed_key_  = false;
    bool found_oid = false;
    while(!found_oid && id_ < mid_) {
      if(ary_->get_key_at_id((start_id_ + id_++) & mask_, *key_, &w_, &o_) == array::FILLED) {
        oid_ = key_->get_bits(0, ary_->lsize());
        found_oid = start_id_ <= oid_ && oid_ < end_id_;
      }
    }

    return found_oid;

  }
};

/// STL like iterator on a large hash array.
template<typename array>
class stl_iterator_base :
    public std::iterator<std::forward_iterator_tag, typename array::value_type>,
    public array::lazy_iterator
{
public:
  typedef typename array::key_type    key_type;
  typedef typename array::mapped_type mapped_type;
  typedef typename array::value_type  value_type;

protected:
  typedef typename array::lazy_iterator           lit;
  typedef std::pair<key_type&, mapped_type> pair;
  pair val_;

public:
  explicit stl_iterator_base(const array* ary, size_t start_id = 0) :
    lit(ary, start_id, ary->size), val_(lit::key_, (mapped_type)0)
  { ++*this; }
  stl_iterator_base(const array* ary, size_t start_id, size_t end_id) :
    lit(ary, start_id, end_id), val_(lit::key_, (mapped_type)0)
  { ++*this; }
  explicit stl_iterator_base() : lit(0, 0, 0), val_(lit::key_, (mapped_type)0) { }
  stl_iterator_base(const stl_iterator_base& rhs) : lit(rhs), val_(lit::key_, rhs.val_.second) { }

  bool operator==(const stl_iterator_base& rhs) const { return lit::ary_ == rhs.ary_ && lit::id_ == rhs.id_; }
  bool operator!=(const stl_iterator_base& rhs) const { return !(*this == rhs); }

  const value_type& operator*() {
    lit::key();
    val_.second = lit::val();
    return val_;
  }
  const value_type* operator->() { return &this->operator*(); }

  stl_iterator_base& operator++() {
    if(!lit::next()) {
      lit::ary_ = 0;
      lit::id_  = 0;
    }
    return *this;
  }

  stl_iterator_base operator++(int) {
    stl_iterator_base res(*this);
    ++*this;
    return res;
  }
};
} } // namespace jellyfish { namespace large_hash {
#endif /* __LARGE_HASH_ITERATOR_HPP__ */
