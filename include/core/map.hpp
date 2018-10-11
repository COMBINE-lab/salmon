#ifndef CORE_MAP_HPP
#define CORE_MAP_HPP

#include <core/memory_resource.hpp>
#include <core/functional.hpp>

#include <map>

namespace core {
inline namespace v2 {

template <class K, class V, class C, class A, class Predicate>
void erase_if (::std::multimap<K, T, C, A>& m, Predicate pred) {
  for (auto iter = begin(m); iter != end(m);) {
    invoke(pred, *iter) ? iter = m.erase(iter) : ++iter;
  }
}

template <class K, class V, class C, class A, class Predicate>
void erase_if (::std::map<K, T, C, A>& m, Predicate pred) {
  for (auto iter = begin(m); iter != end(m);) {
    invoke(pred, *iter) ? iter = m.erase(iter) : ++iter;
  }
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class Key, class T, class Compare=::std::less<Key>>
using multimap = ::std::multimap<
  Key,
  T,
  Compare,
  polymorphic_allocator<::std::pair<Key const, T>>
>;

template <class Key, class T, class Compare=::std::less<Key>>
using map = ::std::map<
  Key,
  T,
  Compare,
  polymorphic_allocator<::std::pair<Key const, T>>
>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_MAP_HPP */
