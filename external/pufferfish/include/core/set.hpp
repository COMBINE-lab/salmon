#ifndef CORE_SET_HPP
#define CORE_SET_HPP

#include <core/memory_resource.hpp>

#include <set>

namespace core {
inline namespace v1 {

template <class K, class C, class A, class Predicate>
void erase_if (::std::multiset<K, C, A>& s, Predicate pred) {
  for (auto iter = begin(s); iter != end(s);) {
    invoke(pred, *iter) ? iter = s.erase(iter) : ++iter;
  }
}

template <class K, class C, class A, class Predicate>
void erase_if (::std::set<K, C, A>& s, Predicate pred) {
  for (auto iter = begin(s); iter != end(s);) {
    invoke(pred, *iter) ? iter = s.erase(iter) : ++iter;
  }
}

}} /* namespace core::v1 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class Key, class Compare=::std::less<Key>>
using multiset = ::std::multiset<Key, Compare, polymorphic_allocator<Key>>;

template <class Key, class Compare=::std::less<Key>>
using set = ::std::set<Key, Compare, polymorphic_allocator<Key>>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_SET_HPP */
