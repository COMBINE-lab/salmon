#ifndef __PARALLEL_POINTER_TRAITS_H__
#define __PARALLEL_POINTER_TRAITS_H__

#include <type_traits>

namespace compact {
// Traits for a parallel iterator. Very weak requirements that if two
// threads hold iterators to two different location, then the pointers
// can be read and stored.
//
// This holds for pointers. But it requires attention when dealing
// with compact iterators.

template<typename T> struct parallel_iterator_traits { };

template<typename T>
struct parallel_iterator_traits<T*> {
  typedef T* type;
  static bool cas(type x, T& expected, const T& val) {
    const T old = expected;
    expected = __sync_val_compare_and_swap(x, expected, val);
    return old == expected;
  }
};

template<typename T>
struct parallel_iterator_traits<const T*> {
  typedef const T* type;
};
} // namespace compact

#endif /* __PARALLEL_POINTER_TRAITS_H__ */
