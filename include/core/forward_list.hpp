#ifndef CORE_FORWARD_LIST_HPP
#define CORE_FORWARD_LIST_HPP

#include <core/memory_resource.hpp>
#include <core/functional.hpp>

#include <forward_list>

namespace core {
inline namespace v2 {

template <class T, class A, class Predicate>
void erase_if (::std::forward_list<T, A>& f, Predicate pred) {
  f.remove_if(pred);
}

template <class T, class A, class U>
void erase (::std::forward_list<T, A>& f, U const& value) {
  using ::std::placeholders::_1;
  f.remove_if(::std::bind(equal<>, _1, ::std::cref(value)));
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class T>
using forward_list = ::std::forward_list<T, polymorphic_allocator<T>>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_FORWARD_LIST_HPP */
