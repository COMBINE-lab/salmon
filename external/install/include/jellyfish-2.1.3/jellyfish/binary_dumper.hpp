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

#ifndef __JELLYFISH_BINARY_DUMPER_HPP__
#define __JELLYFISH_BINARY_DUMPER_HPP__

#include <iostream>
#include <cmath>

#include <jellyfish/sorted_dumper.hpp>

namespace jellyfish {
template<typename Key, typename Val>
class binary_writer {
  int val_len_;
  Val max_val_;
  int key_len_;                 // length of output key field in bytes

public:
  binary_writer(int val_len,   // length of value field in bytes
                int key_len) : // length of key field in bits
    val_len_(val_len),
    max_val_(((Val)1 << (8 * val_len)) - 1),
    key_len_(key_len / 8 + (key_len % 8 != 0))
  { }

  int val_len() const { return val_len_; }
  Val max_val() const { return max_val_; }
  int key_len() const { return key_len_; }

  void write(std::ostream& out, const Key& key, const Val val) {
    out.write((const char*)key.data(), key_len_);
    Val v = std::min(max_val_, val);
    out.write((const char*)&v, val_len_);
  }
};

/// Dump a hash array in sorted binary format. The key/value pairs are
/// written in a sorted list according to the hash function order. The
/// k-mer and count are written in binary, byte aligned.
template<typename storage_t>
class binary_dumper : public sorted_dumper<binary_dumper<storage_t>, storage_t> {
  typedef sorted_dumper<binary_dumper<storage_t>, storage_t> super;
  binary_writer<typename super::key_type, uint64_t> writer;

public:
  static const char* format;

  binary_dumper(int val_len, // length of value field in bytes
                int key_len, // length of key field in bits
                int nb_threads, const char* file_prefix,
                file_header* header = 0) :
    super(nb_threads, file_prefix, header),
    writer(val_len, key_len)
  { }

  virtual void _dump(storage_t* ary) {
    if(super::header_) {
      super::header_->update_from_ary(*ary);
      super::header_->format(format);
      super::header_->counter_len(writer.val_len());
    }
    super::_dump(ary);
  }

  void write_key_value_pair(std::ostream& out, typename super::heap_item item) {
    writer.write(out, item->key_, item->val_);
  }
};
template<typename storage_t>
const char* jellyfish::binary_dumper<storage_t>::format = "binary/sorted";

/// Reader of the format written by binary_dumper. Behaves like an
/// iterator (has next() method which behaves similarly to the next()
/// method of the hash array).
/// The header should be of format binary/sorted, but no check is made.
template<typename Key, typename Val>
class binary_reader {
  std::istream&                 is_;
  const int                     val_len_;
  Key                           key_;
  Val                           val_;
  const RectangularBinaryMatrix m_;
  const size_t                  size_mask_;

public:
  binary_reader(std::istream& is, // stream containing data (past any header)
                file_header* header) :  // header which contains counter_len, matrix, size and key_len
    is_(is), val_len_(header->counter_len()), key_(header->key_len() / 2),
    m_(header->matrix()),
    size_mask_(header->size() - 1)
  { }

  const Key& key() const { return key_; }
  const Val& val() const { return val_; }
  size_t pos() const { return m_.times(key_) & size_mask_; }

  bool next() {
    key_.template read<1>(is_);
    val_ = 0;
    is_.read((char*)&val_, val_len_);
    return is_.good();
  }
};

template<typename Key, typename Val>
class binary_query_base {
  const char* const             data_;
  const unsigned int            val_len_; // In bytes
  const unsigned int            key_len_; // In bytes
  const RectangularBinaryMatrix m_;
  const size_t                  mask_;
  const size_t                  record_len_;
  const size_t                  last_id_;
  Key                           first_key_, last_key_;
  mutable Key                   mid_key_;
  uint64_t                      first_pos_, last_pos_;

public:
  // key_len passed in bits
  binary_query_base(const char* data, unsigned int key_len, unsigned int val_len, const RectangularBinaryMatrix& m, size_t mask,
                    size_t size) :
    data_(data),
    val_len_(val_len),
    key_len_(key_len / 8 + (key_len % 8 != 0)),
    m_(m),
    mask_(mask),
    record_len_(val_len + key_len_),
    last_id_(size / record_len_),
    first_key_(key_len / 2),
    last_key_(key_len / 2),
    mid_key_(key_len / 2)
  {
    if(size % record_len_ != 0)
      eraise(std::length_error) << "Size of database (" << size << ") must be a multiple of the length of a record ("
                                << record_len_ << ")";
    key_at(0, first_key_);
    first_pos_ = key_pos(first_key_);
    key_at(last_id_ - 1, last_key_);
    last_pos_  = key_pos(last_key_);
  }

  bool val_id(const Key& key,  Val* res, uint64_t* id) const {
    if(last_id_ == 0) return false;
    uint64_t first     = 0;
    uint64_t last      = last_id_;
    uint64_t first_pos = first_pos_;
    uint64_t last_pos  = last_pos_;
    const uint64_t pos = key_pos(key);
    uint64_t cid       = 0;
    if(key == first_key_) goto found;
    cid = last_id_ - 1;
    if(key == last_key_) goto found;
    if(pos < first_pos_ || pos > last_pos_) return false;

    // First a guided binary search
    for(uint64_t diff = last - first; diff >= 8; diff = last - first) {
      cid = first + lrint(diff * ((double)(pos - first_pos) / (double)(last_pos - first_pos)));
      cid = std::max(first + 1, cid);
      cid = std::min(cid, last - 1);
      key_at(cid, mid_key_);
      if(key == mid_key_) goto found;
      uint64_t mid_pos = key_pos(mid_key_);
      if(mid_pos > pos || (mid_pos == pos && mid_key_ > key)) {
        last     = cid;
        last_pos = mid_pos;
      } else {
        first     = cid;
        first_pos = mid_pos;
      }
    }

    // Then a linear search (avoids matrix computation)
    for(cid = first + 1; cid < last; ++cid) {
      key_at(cid, mid_key_);
      if(key == mid_key_) goto found;
    }
    return false;

  found:
    val_at(cid, res);
    *id = cid;
    return true;
  }

  Val operator[](const Key& key) const {
    Val res;
    uint64_t id;
    if(!val_id(key, &res, &id))
      return 0;
    return res;
  }

  inline Val check(const Key& key) const { return (*this)[key]; }

protected:
  void key_at(size_t id, Key& key) const {
    memcpy(key.data__(), data_ + id * record_len_, key_len_);
    key.clean_msw();
  }
  void val_at(size_t id, Val* val) const {
    *val = 0;
    memcpy(val, data_ + id * record_len_ + key_len_, val_len_);
  }
  uint64_t key_pos(const Key& key) const {
    return m_.times(key) & mask_;
  }
};
}

#endif /* __JELLYFISH_BINARY_DUMPER_HPP__ */
