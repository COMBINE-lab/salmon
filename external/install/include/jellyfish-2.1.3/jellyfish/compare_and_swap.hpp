/*  This file is part of Jflib.

    Jflib is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jflib is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jflib.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _JFLIB_COMPARE_AND_SWAP_H_
#define _JFLIB_COMPARE_AND_SWAP_H_

#include <sys/types.h>
#include <stdint.h>

namespace jflib {
  // // Atomic load (get) and store (set). For now assume that the
  // // architecture does this based on the virtual key word (true for
  // // x86_64). TODO: improve on other architectures.
  // template<typename T>
  // T a_get(const T &x) { return *(volatile T*)&x; }
  // template<typename T, typename U>
  // T &a_set(T &lhs, const U &rhs) {
  //   *(volatile T*)&lhs = rhs;
  //   return lhs;
  // }

  // Numeric type of length rounded up to the size of a
  // word. Undefined, and raise a compilation error, if the length is
  // not a machine word size
  template<typename T, size_t n> union word_t;
  template<typename T> union word_t<T, 1> { typedef uint8_t  w_t; T v; w_t w; };
  template<typename T> union word_t<T, 2> { typedef uint16_t w_t; T v; w_t w; };
  template<typename T> union word_t<T, 4> { typedef uint32_t w_t; T v; w_t w; };
  template<typename T> union word_t<T, 8> { typedef uint64_t w_t; T v; w_t w; };

  /** Type safe version of CAS.
   * @param [in] ptr Memory location.
   * @param [in] ov  Presumed value at location.
   * @param [in] nv  Value to write.
   * @param [out] cv Value at location at time of call.
   * @return true if CAS is successful.
   *
   * The CAS operation is successful if, at the time of call, ov is
   * equal to *ptr, the value at the memory location. In that case, nv
   * is written to *ptr, and when the call returns, cv == ov.
   *
   * If it fails, cv contains *ptr at the time of call.
   */
  template<typename T>
  bool cas(T *ptr, const T &ov, const T &nv, T *cv) {
    typedef word_t<T, sizeof(T)> val_t;
    val_t _cv, _ov, _nv;
    _ov.v = ov;
    _nv.v = nv;
    _cv.w = __sync_val_compare_and_swap((typename val_t::w_t *)ptr, _ov.w, _nv.w);
    *cv = _cv.v;
    return _cv.w == _ov.w;
  }

  /** Type safe version of CAS.  Identical to 4 argument version,
   * except does not return the previous value.
   */
  template<typename T>
  bool cas(T *ptr, const T &ov, const T &nv) {
    T cv;
    return cas(ptr, ov, nv, &cv);
  }
}


#endif /* _JFLIB_COMPARE_AND_SWAP_H_ */
