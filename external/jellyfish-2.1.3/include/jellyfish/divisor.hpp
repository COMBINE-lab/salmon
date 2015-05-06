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



#ifndef __JELLYFISH_DIVISOR_HPP__
#define __JELLYFISH_DIVISOR_HPP__

#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

namespace jflib {
class divisor64 {
  const uint64_t d_;
#ifdef HAVE_INT128
  const uint16_t p_;
  const unsigned __int128 m_;
#endif
  template<typename T>
  static T div_ceil(T x, T y) {
    T q = x / y;
    T r = x % y;
    return q + (r > 0);
  }

  template<typename T>
  static uint16_t ceilLog2(T x, uint16_t r = 0, uint16_t i = 0) {
    if(x > 1)
      return ceilLog2(x >> 1, r + 1, i | (x & 1));
    return r + i;
  }

public:
  explicit divisor64(uint64_t d) :
    d_(d)
#ifdef HAVE_INT128
    , p_(ceilLog2(d_)),
    m_((div_ceil((unsigned __int128)1 << (64 + p_), (unsigned __int128)d_)) & (uint64_t)-1)
#endif
  { }

  divisor64() :
    d_(0)
#ifdef HAVE_INT128
    , p_(0), m_(0)
#endif
  { }

  explicit divisor64(const divisor64& rhs) :
    d_(rhs.d_)
#ifdef HAVE_INT128
    , p_(rhs.p_),
    m_(rhs.m_)
#endif
  { }

  inline uint64_t divide(const uint64_t n) const {
#ifdef HAVE_INT128
    switch(m_) {
    case 0:
      return n >> p_;
    default:
      const unsigned __int128 n_ = (unsigned __int128)n;
      return (n_ + ((n_ * m_) >> 64)) >> p_;
    }
#else
    return n / d_;
#endif
  }

  inline uint64_t remainder(uint64_t n) const {
#ifdef HAVE_INT128
    switch(m_) {
    case 0:
      return n & (((uint64_t)1 << p_) - 1);
    default:
      return n - divide(n) * d_;
    }
#else
    return n % d_;
#endif
  }

  // Euclidian division: d.division(n, q, r) sets q <- n / d and r
  // <- n % d. This is faster than doing each independently.
  inline void division(uint64_t n, uint64_t &q, uint64_t &r) const {
#ifdef HAVE_INT128
    switch(m_) {
    case 0:
      q = n >> p_;
      r = n & (((uint64_t)1 << p_) - 1);
      break;
    default:
      q = divide(n);
      r = n - q * d_;
      break;
    }
#else
    q = n / d_;
    r = n % d_;
#endif
  }

  uint64_t d() const { return d_; }
  uint64_t p() const {
#ifdef HAVE_INT128
    return p_;
#else
    return 0;
#endif
  }
  uint64_t m() const {
#ifdef HAVE_INT128
    return m_;
#else
    return 0;
#endif
  }
};

inline uint64_t operator/(uint64_t n, const divisor64& d) {
  return d.divide(n);
}
inline uint64_t operator%(uint64_t n, const divisor64& d) {
  return d.remainder(n);
}

inline std::ostream& operator<<(std::ostream& os, const divisor64& d) {
  return os << "d:" << d.d() << ",p:" << d.p() << ",m:" << d.m();
}

} // namespace jflib

#endif /* __JELLYFISH_DIVISOR_HPP__ */

