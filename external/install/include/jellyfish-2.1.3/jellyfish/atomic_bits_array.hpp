/* Quorum
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __JELLYFISH_ATOMIC_BITS_ARRAY_HPP__
#define __JELLYFISH_ATOMIC_BITS_ARRAY_HPP__

#include <stdexcept>
#include <iterator>

#include <jellyfish/allocators_mmap.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/divisor.hpp>

namespace jellyfish {
template<typename Value, typename T, typename Derived>
class atomic_bits_array_base {
  static const int       w_ = sizeof(T) * 8;
  const int              bits_;
  const size_t           size_;
  const T                mask_;
  const jflib::divisor64 d_;
  size_t                 size_bytes_;
  T*                     data_;
  static atomic::gcc     atomic_;

  friend class iterator;
  class iterator : public std::iterator<std::input_iterator_tag, Value> {
    friend class atomic_bits_array_base;
    const atomic_bits_array_base& ary_;
    T*                            word_;
    T                             mask_;
    int                           off_;

    iterator(const atomic_bits_array_base& a, T* w, T m, int o) : ary_(a), word_(w), mask_(m), off_(o) { }
  public:
    bool operator==(const iterator& rhs) const { return word_ == rhs.word_ && off_ == rhs.off_; }
    bool operator!=(const iterator& rhs) const { return word_ != rhs.word_ || off_ != rhs.off_; }
    Value operator*() const { return static_cast<Value>((*word_ & mask_) >> off_); }
    Value* operator->() const { return 0; }
    iterator& operator++() {
      off_ += ary_.bits_;
      if(off_ + ary_.bits_ < w_) {
        mask_ <<= ary_.bits_;
      } else {
        ++word_;
        mask_ = ary_.mask_;
        off_  = 0;
      }
      return *this;
    }
    iterator operator++(int) {
      iterator res(*this);
      ++*this;
      return res;
    }
  };

  class element_proxy {
    T*        word_;
    const T   mask_;
    const int off_;
    T         prev_word_;

    Value get_val(T v) const {
      return static_cast<Value>((v & mask_) >> off_);
    }

  public:
    element_proxy(T* word, T mask, int off) :
      word_(word), mask_(mask), off_(off)
    { }

    operator Value() const { return get_val(*word_); }
    Value get() {
      prev_word_ = *word_;
      return get_val(prev_word_);
    }

    bool set(Value& nval) {
      const T new_word    = (prev_word_ & ~mask_) | ((static_cast<T>(nval) << off_) & mask_);
      const T actual_word = atomic_.cas(word_, prev_word_, new_word);
      if(__builtin_expect(actual_word == prev_word_, 1))
        return true;
      prev_word_ = actual_word;
      nval       = get_val(prev_word_);
      return false;
    }
  };

public:
  atomic_bits_array_base(int bits, // Number of bits per entry
                         size_t size) : // Number of entries
    bits_(bits),
    size_(size),
    mask_((T)-1 >> (w_ - bits)), // mask of one entry at the LSB of a word
    d_(w_ / bits),              // divisor of the number of entries per word
    size_bytes_((size / d_ + (size % d_ != 0)) * sizeof(T)),
    data_(static_cast<Derived*>(this)->alloc_data(size_bytes_))
  {
    static_assert(sizeof(T) >= sizeof(Value), "Container type T must have at least as many bits as value type");
    if((size_t)bits > sizeof(Value) * 8)
      throw std::runtime_error("The number of bits per entry must be less than the number of bits in the value type");
    if(!data_)
      throw std::runtime_error("Can't allocate memory for atomic_bits_array");
  }

  // Return the element at position pos. No check for out of bounds.
  element_proxy operator[](size_t pos) {
    uint64_t q, r;
    d_.division(pos, q, r);
    const int off = r * bits_;
    return element_proxy(data_ + q, mask_ << off, off);
  }
  const element_proxy operator[](size_t pos) const {
    uint64_t q, r;
    d_.division(pos, q, r);
    const int off = r * bits_;
    return element_proxy(data_ + q, mask_ << off, off);
  }
  void write(std::ostream& os) const {
    os.write((const char*)data_, size_bytes_);
  }
  size_t size_bytes() const { return size_bytes_; }
  int bits() const { return bits_; }

  iterator begin() const { return iterator(*this, data_, mask_, 0); }
  iterator end() const {
    uint64_t q, r;
    d_.division(size_, q, r);
    const int off = r * bits_;
    return iterator(*this, data_ + q, mask_ << off, off);
  }
};

template<typename Value, typename T = uint64_t>
class atomic_bits_array :
    protected allocators::mmap,
    public atomic_bits_array_base<Value, T, atomic_bits_array<Value, T> >
{
  typedef atomic_bits_array_base<Value, T, atomic_bits_array<Value, T> > super;
  friend class atomic_bits_array_base<Value, T, atomic_bits_array<Value, T> >;
public:
  atomic_bits_array(int bits, size_t size) :
    allocators::mmap(),
    super(bits, size)
  { }

protected:
  T* alloc_data(size_t s) {
    allocators::mmap::realloc(s);
    return (T*)allocators::mmap::get_ptr();
  }
};

struct mem_info {
  void*  ptr_;
  size_t bytes_;
  mem_info(void* ptr, size_t bytes) : ptr_(ptr), bytes_(bytes) { }
};
template<typename Value, typename T = uint64_t>
class atomic_bits_array_raw :
    protected mem_info,
    public atomic_bits_array_base<Value, T, atomic_bits_array<Value, T> >
{
  typedef atomic_bits_array_base<Value, T, atomic_bits_array<Value, T> > super;
  friend class atomic_bits_array_base<Value, T, atomic_bits_array<Value, T> >;
public:
  atomic_bits_array_raw(void* ptr, size_t bytes, int bits, size_t size) :
    mem_info(ptr, bytes),
    super(bits, size)
  { }

protected:
  T* alloc_data(size_t s) {
    assert(bytes_ == s);
    return (T*)ptr_;
  }
};

} // namespace jellyfish

#endif /* __JELLYFISH_ATOMIC_BITS_ARRAY_HPP__ */
