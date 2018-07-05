#ifndef CORE_DEQUE_HPP
#define CORE_DEQUE_HPP

#include <core/memory_resource.hpp>
#include <core/algorithm.hpp>

#include <deque>

namespace core {
inline namespace v2 {

template <class T, class A, class Predicate>
void erase_if (::std::deque<T, A>& deq, Predicate pred) {
  deq.erase(::core::remove_if(deq, pred), end(deq));
}

template <class T, class A, class U>
void erase (::std::deque<T, A>& deq, U const& value) {
  deq.erase(::core::remove(deq, value), end(deq));
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class T> using deque = ::std::deque<T, polymorphic_allocator<T>>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_DEQUE_HPP */
