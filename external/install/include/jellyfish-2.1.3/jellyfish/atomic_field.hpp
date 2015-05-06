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


#ifndef __JELLYFISH_ATOMIC_FIELD_HPP__
#define __JELLYFISH_ATOMIC_FIELD_HPP__

#include <jellyfish/compare_and_swap.hpp>

namespace jflib {
  /* Define a_get, a_set and a_update
   */
  template <typename T>
  T a_load(T *x) { return *(volatile T*)x; }
  template<typename T, typename U>
  T a_store(T* lhs, const U& rhs) {
    return (*(volatile T*)lhs = rhs);
  }
  template<typename T>
  T* a_load_ptr(T* x) { return a_load((T**)&x); }
  template<typename T, typename U>
  T* a_store_ptr(T* x, const U& rhs) { return a_store((T**)&x, rhs); }

  /** Set value to f(value).
   * @return f(value)
   *
   * The function f may be called more than once. Be careful about
   * side effects (probably better if f has no side effects).
   */
  template<typename T>
  T a_update(T* x, T (*f)(T)) {
    T ov(a_load(x));
    T nv(f(ov));
    while(!cas(x, ov, nv, &ov)) { nv = f(ov); }
    return nv;
  }
  template<typename T>
  T a_load(T &x) { return a_load(&x); }
  template<typename T, typename U>
  T a_store(T &lhs, const U& rhs) { return a_store(&lhs, rhs); }

  /* POD with atomic operators.
   */
  template<typename T>
  struct atomic_pod {
    typedef T type;
    T x;
  };

#define AF_COMPOUND_ASSIGN(op)                                          \
  template<typename T, typename U>                                      \
  T operator op ## = (atomic_pod<T> &x, const U &rhs) { \
    T ov(a_load(&x.x));                                 \
    T nv(ov op rhs);                                                    \
    while(!cas(&x.x, ov, nv, &ov)) { nv = ov op rhs; }                  \
    return nv;                                                          \
  }
  AF_COMPOUND_ASSIGN(+);
  AF_COMPOUND_ASSIGN(-);
  AF_COMPOUND_ASSIGN(*);
  AF_COMPOUND_ASSIGN(/);
  AF_COMPOUND_ASSIGN(%);
  AF_COMPOUND_ASSIGN(>>);
  AF_COMPOUND_ASSIGN(<<);
  AF_COMPOUND_ASSIGN(&);
  AF_COMPOUND_ASSIGN(|);
  AF_COMPOUND_ASSIGN(^);

  /** Set value to f(value).
   * @return f(value)
   *
   * The function f may be called more than once. Be careful about
   * side effects (probably better if f has no side effects).
   */
  template<typename T>
  T a_load(atomic_pod<T> &x) { return a_load(&x.x); }
  template<typename T, typename U>
  T a_store(atomic_pod<T> &lhs, const U &rhs) {
    return a_store(&lhs.x, rhs);
  }
  template<typename T>
  T a_update(atomic_pod<T> &x, T (*f)(T)) {
    return a_update(&x.x, f);
  }

  /* Similar to an atomic_pod, but not a POD, because of its
     constructor and other member functions. Easier to use.
   */
  template<typename T>
  class atomic_field : public atomic_pod<T> {
  public:
    typedef typename atomic_pod<T>::type type;
    explicit atomic_field() { }
    explicit atomic_field(const T& v) { a_store(&this->x, v); }
    atomic_field& operator=(const atomic_pod<T>& rhs) { a_store(&this->x, rhs.x); return *this; }
    atomic_field& operator=(const T& v) { a_store(&this->x, v); return *this; }
    operator T() const { return a_load(&this->x); }
    T update(T (*f)(T)) { return a_update(&this->x, f); }
  };

  template<typename T>
  T a_load(atomic_field<T> &x) { return a_load((atomic_pod<T>&)x); }
  template<typename T, typename U>
  atomic_field<T>& a_store(atomic_field<T>& lhs, const U& rhs) { a_store((atomic_pod<T>&)lhs, rhs); return lhs; }
  template<typename T>
  T a_update(atomic_field<T>& x, T (*f)(T)) { return a_update((atomic_pod<T>&)x, f); }


  /* Allows atomic operation on any (already allocated) data.
   */
  template<typename T>
  class atomic_ref {
    T *ptr;
  public:
    typedef T type;
    explicit atomic_ref(T& x) : ptr(&x) { }
    explicit atomic_ref(T* x) : ptr(x) { }
    atomic_ref& operator=(const T& v) { a_store(ptr, v); return *this; }
    operator T() const { assert(ptr != 0); return a_load(ptr); }
    T* operator&() const { return ptr; }
  };

#define AR_COMPOUND_ASSIGN(op)                          \
  template<typename T, typename U>                      \
  T operator op ## = (atomic_ref<T> &x, const U &rhs) { \
    T ov(x);                                            \
    T nv(ov op rhs);                                    \
    while(!cas(&x, ov, nv, &ov)) { nv = ov op rhs; }    \
    return nv;                                          \
  }
  AR_COMPOUND_ASSIGN(+);
  AR_COMPOUND_ASSIGN(-);
  AR_COMPOUND_ASSIGN(*);
  AR_COMPOUND_ASSIGN(/);
  AR_COMPOUND_ASSIGN(%);
  AR_COMPOUND_ASSIGN(>>);
  AR_COMPOUND_ASSIGN(<<);
  AR_COMPOUND_ASSIGN(&);
  AR_COMPOUND_ASSIGN(|);
  AR_COMPOUND_ASSIGN(^);
}


#endif /* __JELLYFISH_ATOMIC_FIELD_HPP__ */
