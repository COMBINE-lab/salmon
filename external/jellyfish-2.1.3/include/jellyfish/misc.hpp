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

#ifndef __JELLYFISH_MISC_HPP__
#define __JELLYFISH_MISC_HPP__

#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <exception>
#include <stdexcept>
#include <string>
#include <new>
#include <ostream>
#include <utility>
#include <iterator>
#include <algorithm>

namespace jellyfish {
#define bsizeof(v)      (8 * sizeof(v))
typedef uint_fast64_t uint_t;
//#define UINT_C(x)
#define PRIUINTu PRIuFAST64
#define PRIUINTx PRIxFAST64

inline int leading_zeroes(int x) { return __builtin_clz(x); } // CLK
inline int leading_zeroes(unsigned int x) { return __builtin_clz(x); }
inline int leading_zeroes(unsigned long x) { return __builtin_clzl(x); }
inline int leading_zeroes(unsigned long long x) { return __builtin_clzll(x); }

/// The floor of the log base two of n. Undefined if n == 0
template <typename T>
uint16_t floorLog2(T n) {
  return sizeof(T) * 8 - 1 - leading_zeroes(n);
}

/// The ceiling of the log base two of n. Undefined if n == 0
template<typename T>
uint16_t ceilLog2(T n) {
  uint16_t r = floorLog2(n);
  return n > (((T)1) << r) ? r + 1 : r;
}

/// The ceiling of the quotient of the division of a by b. I.e. if b
/// divides a, then div_ceil(a, b) == a / b. Otherwise, div_ceil(a, b)
/// == a / b + 1
template<typename T>
T div_ceil(T a, T b) {
  T q = a / b;
  return a % b == 0 ? q : q + 1;
}

/// Number of bits necessary to encode number n. Undefined if n ==
/// 0. The following should be true: 2^bitsize(n) - 1 >= n >
/// 2^(bitsize(n) - 1)
template<typename T>
uint16_t bitsize(T n) {
  return floorLog2(n) + 1;
}

inline uint32_t reverse_bits(uint32_t v) {
  // swap odd and even bits
  v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);
  // swap consecutive pairs
  v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);
  // swap nibbles ...
  v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);
  // swap bytes
  v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);
  // swap 2-byte long pairs
  v = ( v >> 16             ) | ( v               << 16);
  return v;
}

inline uint64_t reverse_bits(uint64_t v) {
  v = ((v >> 1)  & 0x5555555555555555UL) | ((v & 0x5555555555555555UL) << 1);
  v = ((v >> 2)  & 0x3333333333333333UL) | ((v & 0x3333333333333333UL) << 2);
  v = ((v >> 4)  & 0x0F0F0F0F0F0F0F0FUL) | ((v & 0x0F0F0F0F0F0F0F0FUL) << 4);
  v = ((v >> 8)  & 0x00FF00FF00FF00FFUL) | ((v & 0x00FF00FF00FF00FFUL) << 8);
  v = ((v >> 16) & 0x0000FFFF0000FFFFUL) | ((v & 0x0000FFFF0000FFFFUL) << 16);
  v = ( v >> 32                        ) | ( v                         << 32);
  return v;
}

uint64_t bogus_sum(void *data, size_t len);

template <typename T>
size_t bits_to_bytes(T bits) {
  return (size_t)((bits / 8) + (bits % 8 != 0));
}

template <typename T>
union Tptr {
  void *v;
  T    *t;
};
template <typename T>
T *calloc_align(size_t nmemb, size_t alignment) {
  Tptr<T> ptr;
  if(posix_memalign(&ptr.v, alignment, sizeof(T) * nmemb) < 0)
    throw std::bad_alloc();
  return ptr.t;
}

/* Be pedantic about memory access. Any misaligned access will
 * generate a BUS error.
 */
void disabled_misaligned_mem_access();

/* Raison d'etre of this version of mem_copy: It seems we have slow
 * down due to misaligned cache accesses. glibc memcpy does unaligned
 * memory accesses and crashes when they are disabled. This version
 * does only aligned memory access (see above).
 */
template <typename T>
void mem_copy(char *dest,  const char *src, const T &len) {
  // dumb copying char by char
  for(T i = (T)0; i < len; ++i)
    *dest++ = *src++;
}

/* Slice a large number (total) in almost equal parts. return [start,
   end) corresponding to the ith part (0 <= i < number_of_slices)
 */
template<typename T>
std::pair<T,T> slice(T i, T number_of_slices, T total) {
  const T slice_size   = total / number_of_slices;
  const T slice_remain = total % number_of_slices;

  const T start = std::max((T)0,
                           std::min(total, i * slice_size + (i > 0 ? slice_remain : 0)));
  const T end   = std::max((T)0,
                           std::min(total, (i + 1) * slice_size + slice_remain));

  return std::make_pair(start, end);
}

uint64_t random_bits(int length);
inline uint64_t random_bits() { return random_bits(64); }

// Quote string that could contain shell special characters
std::string quote_arg(const std::string& arg);

std::streamoff get_file_size(std::istream& is);

/// Find the first element for which the predicate p is false. The
/// input range [first, last) is assumed to be sorted according to the
/// predicate p: p(x) is false and p(y) is true implies x comes after
/// y in the input range. (I.e., the elements for which p(x) is true
/// come first followed by the elements for which p(x) is false).
template <class ForwardIterator, class Predicate>
ForwardIterator binary_search_first_false(ForwardIterator first, ForwardIterator last, Predicate p)
{
  ForwardIterator it;
  typename std::iterator_traits<ForwardIterator>::difference_type count, step;
  count = std::distance(first,last);
  while(count>0) {
    it = first; step = count / 2; std::advance(it,step);
    if(p(*it)) {
      first  = ++it;
      count -= step + 1;
    } else
      count=step;
  }
  return first;
}

/// An integer type which behaves like a random pointer to
/// itself. Meaning, with `it(5)`, then `*it == 5` and `*++it ==
/// 6`. In other words, it is a pointer to an array `a` initialized
/// with `a[i] = i`, except the array is not instantiated and does not
/// have a fixed size.
template<typename T>
class pointer_integer : public std::iterator<std::random_access_iterator_tag, T> {
    T x_;
  typedef typename std::iterator<std::random_access_iterator_tag, T> super;
 public:
  typedef T                                 value_type;
  typedef typename super::difference_type   difference_type;
  typedef typename super::pointer           pointer;
  typedef typename super::reference         reference;
  typedef typename super::iterator_category iterator_category;

  pointer_integer() : x_(0) { }
  explicit pointer_integer(T x) : x_(x) { }
  pointer_integer(const pointer_integer& rhs) : x_(rhs.x_) { }
  pointer_integer& operator=(const pointer_integer& rhs) {
    x_ = rhs.x_;
    return *this;
  }
  pointer_integer& operator++() { ++x_; return *this; }
  pointer_integer operator++(int) {
    pointer_integer res(*this);
    ++x_;
    return res;
  }
  pointer_integer& operator--() { --x_; return *this; }
  pointer_integer operator--(int) {
    pointer_integer res(*this);
    --x_;
    return res;
  }

  bool operator==(const pointer_integer& rhs) const { return x_ == rhs.x_; }
  bool operator!=(const pointer_integer& rhs) const { return x_ != rhs.x_; }
  bool operator<(const pointer_integer& rhs) const { return x_ < rhs.x_; }
  bool operator>(const pointer_integer& rhs) const { return x_ > rhs.x_; }
  bool operator<=(const pointer_integer& rhs) const { return x_ <= rhs.x_; }
  bool operator>=(const pointer_integer& rhs) const { return x_ >= rhs.x_; }

  reference operator*() { return x_; }
  pointer operator->() { return &x_; } // Probably useless

  difference_type operator-(pointer_integer& rhs) { return x_ - rhs.x_; }

  pointer_integer operator+(T x) const { return pointer_integer(x_ + x); }
  pointer_integer operator-(T x) const { return pointer_integer(x_ - x); }
  pointer_integer& operator+=(T x) { x_ += x; return *this; }
  pointer_integer& operator-=(T x) { x_ -= x; return *this; }

  value_type operator[](T i) const { return x_ + i; }
};

template<typename T>
pointer_integer<T> operator+(T x, pointer_integer<T>& p) { return pointer_integer<T>(x + *p); }
template<typename T>
pointer_integer<T> operator-(T x, pointer_integer<T>& p) { return pointer_integer<T>(x - *p); }
} // namespace jellyfish

#endif // __MISC_HPP__
