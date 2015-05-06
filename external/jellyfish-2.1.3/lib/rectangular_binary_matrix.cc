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

#include <new>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <assert.h>
#include <jellyfish/rectangular_binary_matrix.hpp>

uint64_t *jellyfish::RectangularBinaryMatrix::alloc(unsigned int r, unsigned int c) {
  if(r > (sizeof(uint64_t) * 8) || r == 0 || c == 0) {
    std::ostringstream err;
    err << "Invalid matrix size " << r << "x" << c;
    throw std::out_of_range(err.str());
  }
  void *mem;
  // Make sure the number of words allocated is a multiple of
  // 8. Necessary for loop unrolling of vector multiplication
  size_t alloc_columns = (c / 8 + (c % 8 != 0)) * 8;
  if(posix_memalign(&mem, sizeof(uint64_t) * 2, alloc_columns * sizeof(uint64_t)))
    throw std::bad_alloc();
  memset(mem, '\0', sizeof(uint64_t) * alloc_columns);
  return (uint64_t *)mem;
}

void jellyfish::RectangularBinaryMatrix::init_low_identity() {
  memset(_columns, '\0', sizeof(uint64_t) * _c);
  unsigned int row = std::min(_c, _r);
  unsigned int col = _c - row;
  _columns[col] = (uint64_t)1 << (row - 1);
  for(unsigned int i = col + 1; i < _c; ++i)
    _columns[i] = _columns[i - 1] >> 1;
}

bool jellyfish::RectangularBinaryMatrix::is_low_identity() {
  unsigned int row = std::min(_c, _r);
  unsigned int col = _c - row;

  for(unsigned int i = 0; i < col; ++i)
    if(_columns[i])
      return false;
  if(_columns[col] != (uint64_t)1 << (row - 1))
    return false;
  for(unsigned int i = col + 1; i < _c; ++i)
    if(_columns[i] != _columns[i - 1] >> 1)
      return false;
  return true;
}

jellyfish::RectangularBinaryMatrix jellyfish::RectangularBinaryMatrix::pseudo_multiplication(const jellyfish::RectangularBinaryMatrix &rhs) const {
  if(_r != rhs._r || _c != rhs._c)
    throw std::domain_error("Matrices of different size");
  RectangularBinaryMatrix res(_r, _c);

  // v is a vector. The lower part is equal to the given column of rhs
  // and the high part is the identity matrix.
  //  uint64_t v[nb_words()];
  uint64_t *v = new uint64_t[nb_words()];
  memset(v, '\0', sizeof(uint64_t) * nb_words());
  unsigned int j = nb_words() - 1;
  v[j]           = msb();
  const unsigned int row = std::min(_c, _r);
  const unsigned int col = _c - row;

  unsigned int i;
  for(i = 0; i < col; ++i) {
    // Set the lower part to rhs and do vector multiplication
    v[0] ^= rhs[i];
    res.get(i) = this->times(&v[0]);
    //res.get(i) = this->times_loop(v);

    // Zero the lower part and shift the one down the diagonal.
    v[0]  ^= rhs[i];
    v[j] >>= 1;
    if(!v[j])
      v[--j] = (uint64_t)1 << (sizeof(uint64_t) * 8 - 1);
  }
  // No more identity part to deal with
  memset(v, '\0', sizeof(uint64_t) * nb_words());
  for( ; i < _c; ++i) {
    v[0] = rhs[i];
    res.get(i) = this->times(v);
    //res.get(i) = this->times_loop(v);
  }

  delete[] v;
  return res;
}

unsigned int jellyfish::RectangularBinaryMatrix::pseudo_rank() const {
  unsigned int            rank = _c;
  RectangularBinaryMatrix pivot(*this);

  // Make the matrix lower triangular.
  unsigned int srow = std::min(_r, _c);
  unsigned int scol = _c - srow;
  uint64_t mask = (uint64_t)1 << (srow - 1);
  for(unsigned int i = scol; i < _c; ++i, mask >>= 1) {
    if(!(pivot.get(i) & mask)) {
      // current column has a 0 in the diagonal. XOR it with another
      // column to get a 1.
      unsigned int j;
      for(j = i + 1; j < _c; ++j)
        if(pivot.get(j) & mask)
          break;
      if(j == _c) {
        // Did not find one, the matrix is not full rank.
        rank = i;
        break;
      }
      pivot.get(i) ^= pivot.get(j);
    }

    // Zero out every ones on the ith row in the upper part of the
    // matrix.
    for(unsigned int j = i + 1; j < _c; ++j)
      if(pivot.get(j) & mask)
        pivot.get(j) ^= pivot.get(i);
  }

  return rank;
}

jellyfish::RectangularBinaryMatrix jellyfish::RectangularBinaryMatrix::pseudo_inverse() const {
  RectangularBinaryMatrix pivot(*this);
  RectangularBinaryMatrix res(_r, _c); res.init_low_identity();
  unsigned int            i, j;
  uint64_t                mask;

  // Do gaussian elimination on the columns and apply the same
  // operation to res.

  // Make pivot lower triangular.
  unsigned int srow = std::min(_r, _c);
  unsigned int scol = _c - srow;
  mask = (uint64_t)1 << (srow - 1);
  for(i = scol; i < _c; ++i, mask >>= 1) {
    if(!(pivot.get(i) & mask)) {
      // current column has a 0 in the diagonal. XOR it with another
      // column to get a 1.
      unsigned int j;
      for(j = i + 1; j < _c; ++j)
        if(pivot.get(j) & mask)
          break;
      if(j == _c)
        throw std::domain_error("Matrix is singular");
      pivot.get(i) ^= pivot.get(j);
      res.get(i)   ^= res.get(j);
    }
    // Zero out every ones on the ith row in the upper part of the
    // matrix.
    for(j = i + 1; j < _c; ++j) {
      if(pivot.get(j) & mask) {
        pivot.get(j) ^= pivot.get(i);
        res.get(j)   ^= res.get(i);
      }
    }
  }

  // Make pivot the lower identity
  mask = (uint64_t)1 << (srow - 1);
  for(i = scol; i < _c; ++i, mask >>= 1) {
    for(j = 0; j < i; ++j) {
      if(pivot.get(j) & mask) {
        pivot.get(j) ^= pivot.get(i);
        res.get(j)   ^= res.get(i);
      }
    }
  }

  return res;
}

void jellyfish::RectangularBinaryMatrix::print(std::ostream &os) const {
  uint64_t mask = (uint64_t)1 << (_r - 1);
  for( ; mask; mask >>= 1) {
    for(unsigned int j = 0; j < _c; ++j) {
      os << (mask & _columns[j] ? "1" : "0");
    }
    os << "\n";
  }
}

template<typename T>
void jellyfish::RectangularBinaryMatrix::print_vector(std::ostream &os, const T &v) const {
  uint64_t mask = msb();
  for(int i = nb_words() - 1; i >= 0; --i) {
    for( ; mask; mask >>= 1)
      os << (v[i] & mask ? "1" : "0");
    mask  = (uint64_t)1 << (sizeof(uint64_t) * 8 - 1);
  }
  os << "\n";
}

jellyfish::RectangularBinaryMatrix jellyfish::RectangularBinaryMatrix::randomize_pseudo_inverse(uint64_t (*rng)()) {
  while(true) {
    randomize(rng);
    try {
      return pseudo_inverse();
    } catch(std::domain_error &e) { }
  }
}
