#ifndef __CONST_ITERATOR_TRAITS_H__
#define __CONST_ITERATOR_TRAITS_H__

#include <type_traits>


namespace compact {
template<typename T> struct const_iterator_traits { };

template<typename T>
struct const_iterator_traits<T*> {
  typedef typename std::add_const<T>::type* type;
};

template<typename T>
struct const_iterator_traits<const T*> {
  typedef const T* type;
};
} // namespace compact

#endif /* __CONST_ITERATOR_TRAITS_H__ */
