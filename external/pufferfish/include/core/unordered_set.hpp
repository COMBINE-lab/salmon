#ifndef CORE_UNORDERED_SET_HPP
#define CORE_UNORDERED_SET_HPP

#include <core/memory_resource.hpp>

#include <unordered_set>

namespace core {
inline namespace v2 {

template <class K, class H, class P, class A, class Predicate>
void erase_if (::std::unordered_multiset<K, H, P, A>& m, Predicate pred) {
  for (auto iter = begin(m); iter != end(m);) {
    invoke(pred, *iter) ? iter = m.erase(iter) : ++iter;
  }
}

template <class K, class H, class P, class A, class Predicate>
void erase_if (::std::unordered_set<K, H, P, A>& m, Predicate pred) {
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
  class Hash=::std::hash<Key>,
  class Pred=::std::equal_to<Key>
> using unordered_multiset = ::std::unordered_multiset<
  Key,
  Hash,
  Pred,
  polymorphic_allocator<Key>
>;

template <
  class Key,
  class Hash=::std::hash<Key>,
  class Pred=::std::equal_to<Key>
> using unordered_set = ::std::unordered_set<
  Key,
  Hash,
  Pred,
  polymorphic_allocator<Key>>
>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_UNORDERED_SET_HPP */
