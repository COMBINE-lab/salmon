#ifndef __PREFETCH_TRAITS_H__
#define __PREFETCH_TRAITS_H__

#include <type_traits>

namespace compact {
// Traits to prefetch an iterator

template<typename T> struct prefetch_iterator_traits { };

template<typename T>
struct prefetch_iterator_traits<T*> {
  template<int level = 0>
  static void read(T* ptr) { __builtin_prefetch((void*)ptr, 0, level); }
  template<int level = 0>
  static void write(T* ptr) { __builtin_prefetch((void*)ptr, 1, level); }
};

template<typename T>
struct prefetch_iterator_traits<const T*> {
  template<int level = 0>
  static void read(const T* ptr) { __builtin_prefetch((void*)ptr, 0, level); }
  template<int level = 0>
  static void write(const T* ptr) { __builtin_prefetch((void*)ptr, 1, level); }
};
} // namespace compact

#endif /* __PREFETCH_TRAITS_H__ */
