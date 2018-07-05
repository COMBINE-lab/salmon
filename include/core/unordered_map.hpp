#ifndef CORE_UNORDERED_MAP_HPP
#define CORE_UNORDERED_MAP_HPP

#include <core/memory_resource.hpp>

#include <unordered_map>

namespace core {
inline namespace v2 {

template <class K, class T, class H, class P, class A, class Predicate>
void erase_if (::std::unordered_multimap<K, T, H, P, A>& m, Predicate pred) {
  for (auto iter = begin(m); iter != end(m);) {
    invoke(pred, *iter) ? iter = m.erase(iter) : ++iter;
  }
}

template <class K, class T, class H, class P, class A, class Predicate>
void erase_if (::std::unordered_map<K, T, H, P, A>& m, Predicate pred) {
  for (auto iter = begin(m); iter != end(m);) {
    invoke(pred, *iter) ? iter = m.erase(iter) : ++iter;
  }
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <
  class Key,
  class T,
  class Hash=::std::hash<Key>,
  class Pred=::std::equal_to<Key>
> using unordered_multimap = ::std::unordered_multimap<
  Key,
  T,
  Hash,
  Pred,
  polymorphic_allocator<::std::pair<Key const, T>>
>;

template <
  class Key,
  class T,
  class Hash=::std::hash<Key>,
  class Pred=::std::equal_to<Key>
> using unordered_map = ::std::unordered_map<
  Key,
  T,
  Hash,
  Pred,
  polymorphic_allocator<::std::pair<Key const, T>>
>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_UNORDERED_MAP_HPP */
