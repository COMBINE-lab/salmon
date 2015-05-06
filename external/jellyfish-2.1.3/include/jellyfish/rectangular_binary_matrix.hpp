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

#ifndef __JELLYFISH_RECTANGULAR_BINARY_MATRIX_HPP__
#define __JELLYFISH_RECTANGULAR_BINARY_MATRIX_HPP__

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <jellyfish/misc.hpp>
#include <iostream>
#include <exception>
#include <stdexcept>
#include <vector>
#include <limits>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

// Column major representation
//
// Rectangular matrices on Z/2Z of size _r x _c where 1<=_r<=64 (_c is
// not limited) and _r <= _c. I.e., the matrix can be stored as an
// array of 64 bit word, each representing a column (the highest 64-_r
// bits of each word are set to 0).
//
// Multiplication between a matrix and vector of size _c x 1 gives a
// vector of size _r x 1 stored as one 64 bit word.

namespace jellyfish {
  class RectangularBinaryMatrix {
  public:
    RectangularBinaryMatrix(unsigned int r, unsigned c)
      : _columns(alloc(r, c)), _r(r), _c(c) { }
    RectangularBinaryMatrix(const RectangularBinaryMatrix &rhs)
    : _columns(alloc(rhs._r, rhs._c)), _r(rhs._r), _c(rhs._c) {
      memcpy(_columns, rhs._columns, sizeof(uint64_t) * _c);
    }
    RectangularBinaryMatrix(RectangularBinaryMatrix&& rhs) :
    _columns(rhs._columns), _r(rhs._r), _c(rhs._c) {
      rhs._columns = 0;
    }
    // Initialize from raw data. raw must contain at least c words.
    template<typename T>
    RectangularBinaryMatrix(const T &raw, unsigned int r, unsigned c)
      : _columns(alloc(r, c)), _r(r), _c(c) {
      for(unsigned int i = 0; i < _c; ++i)
        _columns[i] = raw[i] & cmask();
    }
    ~RectangularBinaryMatrix() {
      free(_columns);
    }

    RectangularBinaryMatrix &operator=(const RectangularBinaryMatrix &rhs) {
      if(_r != rhs._r || _c != rhs._c)
        throw std::invalid_argument("RHS matrix dimensions do not match");
      memcpy(_columns, rhs._columns, sizeof(uint64_t) * _c);
      return *this;
    }
    RectangularBinaryMatrix& operator=(RectangularBinaryMatrix&& rhs) {
      if(_r != rhs._r || _c != rhs._c)
        throw std::invalid_argument("RHS matrix dimensions do not match");
      std::swap(_columns, rhs._columns);
      return *this;
    }

    bool operator==(const RectangularBinaryMatrix &rhs) const {
      if(_r != rhs._r || _c != rhs._c)
        return false;
      return !memcmp(_columns, rhs._columns, sizeof(uint64_t) * _c);
    }
    bool operator!=(const RectangularBinaryMatrix &rhs) const {
      return !(*this == rhs);
    }

    // Get i-th column. No check on range
    const uint64_t & operator[](unsigned int i) const { return _columns[i]; }

    unsigned int r() const { return _r; }
    unsigned int c() const { return _c; }

    // True if every column is zero
    bool is_zero() const {
      uint64_t *p = _columns;
      while(*p == 0 && p < _columns + _c)
        ++p;
      return (p - _columns) == _c;
    }

    // Randomize the content of the matrix
    void randomize(uint64_t (*rng)()) {
      for(unsigned int i = 0; i < _c; ++i)
        _columns[i] = rng() & cmask();
    }
    //void randomize() { randomize(rng); }

    // Make and check that the matrix the lower right corner of the
    // identity.
    void init_low_identity();
    bool is_low_identity();

    // Left matrix vector multiplication. Type T supports the operator
    // v[i] to return the i-th 64 bit word of v.
    template<typename T>
    uint64_t times_loop(const T &v) const;


#ifdef HAVE_SSE
    // This SSE implementation only works if the number of columns is
    // even.
    template<typename T>
    uint64_t times_sse(const T &v) const;
#endif

#ifdef HAVE_INT128
    // Implementation using __int128
    template<typename T>
    uint64_t times_128(const T& v) const;
#endif

    template<typename T>
    inline uint64_t times(const T& v) const {
#ifdef HAVE_SSE
      return times_sse(v);
#elif HAVE_INT128
      return times_128(v);
#else
      return times_loop(v);
#endif
    }

    // Return a matrix which is the "pseudo inverse" of this matrix. It
    // is assumed that there is above this square matrix an identity
    // block and a zero so as to make the matrix squared. Raise an
    // exception std::domain_error if the matrix is singular.
    RectangularBinaryMatrix pseudo_inverse() const;

    // Return the multiplication of this and rhs. As in pseudo_inverse,
    // the two matrices are viewed as being squared, padded above by the
    // identity.
    RectangularBinaryMatrix pseudo_multiplication(const RectangularBinaryMatrix &rhs) const;

    // Initialize the object with a pseudo-invertible matrix and return its pseudo-inverse
    RectangularBinaryMatrix randomize_pseudo_inverse(uint64_t (*rng)());
    RectangularBinaryMatrix randomize_pseudo_inverse() { return randomize_pseudo_inverse(random_bits); }

    // Return the rank of the matrix. The matrix is assumed to be
    // squared, padded above by the identity.
    unsigned int pseudo_rank() const;

    // Display matrix
    void print(std::ostream &os) const;
    template<typename T>
    void print_vector(std::ostream &os, const T &v) const;

    // Nb words in vector for multiplication
    uint64_t nb_words() const { return (_c >> 6) + ((_c & 0x3f) != 0); }
    // Mask of most significant bit in most significant word of a vector
    // with _c rows.
    uint64_t msb() const {
      int shift = _c & 0x3f;
      if(shift == 0)
        shift = sizeof(uint64_t) * 8;
      return (uint64_t)1 << (shift - 1);
    }

  private:
    // Store column by column. A column may use one word.  By
    // convention, the "unused" bits (most significant bits) of each
    // column are set to 0.
    uint64_t *         _columns;
    const unsigned int _r, _c;

    static uint64_t *alloc(unsigned int r, unsigned int c) __attribute__((malloc));
    // Mask for column word (zero msb)
    uint64_t cmask() const { return std::numeric_limits<uint64_t>::max() >> (std::numeric_limits<uint64_t>::digits - _r); }
    // Mask of highest word of a vector with _c rows (Most Significant
    // Word)
    uint64_t msw() const { return (msb() << 1) - 1; }
    // Nb of bits used in highest word of vector with _c rows.
    uint64_t nb_msb() const {
      uint64_t nb = _c & 0x3f;
      return nb ? nb : sizeof(uint64_t) * 8;
    }
    // Allow to change the matrix vectors. No check on i.
    uint64_t & get(unsigned int i) { return _columns[i]; }
  };

  template<typename T>
  uint64_t RectangularBinaryMatrix::times_loop(const T &v) const {
    uint64_t       *p   = _columns + _c - 1;
    uint64_t        res = 0, x = 0, j = 0;
    const uint64_t  one = (uint64_t)1;

    for(unsigned int i = 0; i < nb_words(); ++i) {
      j = sizeof(uint64_t) * 8;
      x = v[i];
      if(i == nb_words() - 1) {
        x &= msw();
        j  = nb_msb();
      }
      for( ; j > 7; j -= 8, p -= 8) {
        res ^= (-(x & one)) & p[0];  x >>= 1;
        res ^= (-(x & one)) & p[-1]; x >>= 1;
        res ^= (-(x & one)) & p[-2]; x >>= 1;
        res ^= (-(x & one)) & p[-3]; x >>= 1;
        res ^= (-(x & one)) & p[-4]; x >>= 1;
        res ^= (-(x & one)) & p[-5]; x >>= 1;
        res ^= (-(x & one)) & p[-6]; x >>= 1;
        res ^= (-(x & one)) & p[-7]; x >>= 1;
      }
    }

    // Finish the loop
    switch(j) {
    case 7: res ^= (-(x & one)) & *p--; x >>= 1;
    case 6: res ^= (-(x & one)) & *p--; x >>= 1;
    case 5: res ^= (-(x & one)) & *p--; x >>= 1;
    case 4: res ^= (-(x & one)) & *p--; x >>= 1;
    case 3: res ^= (-(x & one)) & *p--; x >>= 1;
    case 2: res ^= (-(x & one)) & *p--; x >>= 1;
    case 1: res ^= (-(x & one)) & *p;
    }

    return res;
  }

#ifdef HAVE_SSE
  template<typename T>
  uint64_t RectangularBinaryMatrix::times_sse(const T &v) const {
#define FFs ((uint64_t)-1)
    static const uint64_t smear[8] asm("smear") __attribute__ ((aligned(16),used)) =
      {0, 0, 0, FFs, FFs, 0, FFs, FFs};
    typedef uint64_t xmm_t __attribute__((vector_size(16)));

    uint64_t *p = _columns + _c - 8;

    // //#ifdef __ICC
    // register xmm_t acc;
    // register xmm_t load;
    // memset(&acc, '\0', 16);
    // memset(&load, '\0', 16);
    // #else
#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wuninitialized"
#endif
    register xmm_t acc  = acc ^ acc; // Set acc to 0
    register xmm_t load = load ^ load;
#ifdef __clang__
#pragma clang diagnostic pop
#endif
    // #endif

//     // Zero out acc
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wuninitialized"
//     asm("pxor %0,%0\n\t" : "=x"(acc) : "0"(acc));
//     asm("pxor %0,%0\n\t" : "=x"(load) : "0"(load));
// #pragma GCC diagnostic pop

    // i is the lower 2 bits of x, and an index into the smear array. Compute res ^= smear[i] & p[j].
#define AND_XOR(off)                                                    \
    asm("movdqa (%[s],%[i]), %[load]\n\t"                               \
        "pand " off "(%[p]),%[load]\n\t"                                \
        "pxor %[load],%[acc]\n\t"                                       \
        : [acc]"=&x"(acc)                                               \
        : "[acc]"(acc),  [i]"r"(i), [p]"r"(p), [s]"r"(smear), [load]"x"(load))


    uint64_t i, j = 0, x = 0;
    for(unsigned int w = 0; w < nb_words(); ++w) {
      x = v[w];
      j = sizeof(uint64_t) * 8;
      if(w == nb_words() - 1) {
        x &= msw();
        j  = nb_msb();
      }
      for( ; j > 7; j -= 8, p -= 8) {
        i = (x & (uint64_t)0x3) << 4;
        AND_XOR("0x30");
        x >>= 2;
        i = (x & (uint64_t)0x3) << 4;
        AND_XOR("0x20");
        x >>= 2;
        i = (x & (uint64_t)0x3) << 4;
        AND_XOR("0x10");
        x >>= 2;
        i = (x & (uint64_t)0x3) << 4;
        AND_XOR("");
        x >>= 2;
      }
    }

    // Finish loop
    p = _columns;
    switch(j) {
    case 6:
      i = (x & (uint64_t)0x3) << 4;
      AND_XOR("0x20");
      x >>= 2;
    case 4:
      i = (x & (uint64_t)0x3) << 4;
      AND_XOR("0x10");
      x >>= 2;
    case 2:
      i = (x & (uint64_t)0x3) << 4;
      AND_XOR("");
    }

    // Get result out
    uint64_t res1, res2;
    asm("movd %[acc], %[res1]\n\t"
        "psrldq $8, %[acc]\n\t"
        "movd %[acc], %[res2]\n\t"
        : [res1]"=r"(res1), [res2]"=r"(res2)
        : [acc]"x"(acc));
    return res1 ^ res2;
  }
#endif // HAVE_SSE

#ifdef HAVE_INT128
  template<typename T>
  uint64_t RectangularBinaryMatrix::times_128(const T &v) const {
    typedef unsigned __int128 u128;
    static const u128 smear[4] =
      { (u128)0,
        (((u128)1 << 64) - 1) << 64,
        ((u128)1 << 64) - 1,
        (u128)-1
      };\
    u128  res = res ^ res;
    u128* p   = (u128*)(_columns + _c - 2);

    uint64_t j = 0, x = 0;
    for(unsigned int w = 0; w < nb_words(); ++w) {
      x = v[w];
      j = sizeof(uint64_t) * 8;
      if(w == nb_words() - 1) {
        x &= msw();
        j  = nb_msb();
      }
      for( ; j > 7; j -= 8, p -= 4) {
        res ^= smear[x & (uint64_t)0x3] & p[ 0]; x >>= 2;
        res ^= smear[x & (uint64_t)0x3] & p[-1]; x >>= 2;
        res ^= smear[x & (uint64_t)0x3] & p[-2]; x >>= 2;
        res ^= smear[x & (uint64_t)0x3] & p[-3]; x >>= 2;
      }
    }

    switch(j) {
    case 6: res ^= smear[x & (uint64_t)0x3] & *p--; x >>=2;
    case 4: res ^= smear[x & (uint64_t)0x3] & *p--; x >>=2;
    case 2: res ^= smear[x & (uint64_t)0x3] & *p;
    }

    return (res ^ (res >> 64)) & smear[2];
  }
#endif // HAVE_INT128

}

#endif
