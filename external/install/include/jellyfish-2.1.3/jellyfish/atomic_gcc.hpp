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

#ifndef __JELLYFISH_ATOMIC_GCC_HPP__
#define __JELLYFISH_ATOMIC_GCC_HPP__

namespace atomic
{
  class gcc
  {
  public:
    template<typename T>
    static inline T cas(volatile T *ptr, T oval, T nval) {
      return __sync_val_compare_and_swap(ptr, oval, nval);
    }

    template<typename T>
    static inline T set(T *ptr, T nval) {
      return __sync_lock_test_and_set(ptr, nval);
    }

    template<typename T>
    static inline T add_fetch(volatile T *ptr, T x) {
      T ncount = *ptr, count;
      do {
	count = ncount;
	ncount = cas((T *)ptr, count, count + x);
      } while(ncount != count);
      return count + x;
    }

    template<typename T>
    static inline T fetch_add(volatile T *ptr, T x) {
      T ncount = *ptr, count;
      do {
	count = ncount;
	ncount = cas((T *)ptr, count, (T)(count + x));
      } while(ncount != count);
      return count;
    }

    template<typename T>
    static inline T set_to_max(volatile T *ptr, T x) {
      T count = *ptr;
      while(x > count) {
        T ncount = cas(ptr, count, x);
        if(ncount == count)
          return x;
        count = ncount;
      }
      return count;
    }
  };
}
#endif
