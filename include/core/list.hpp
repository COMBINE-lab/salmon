#ifndef CORE_LIST_HPP
#define CORE_LIST_HPP

#include <core/memory_resource.hpp>
#include <core/functional.hpp>

#include <list>

namespace core {
inline namespace v2 {

template <class T, class A, class Predicate>
void erase_if (::std::list<T, A>& l, Predicate pred) { l.remove_if(pred); }

template <class T, class A, class U>
void erase (::std::list<T, A>& l, U const& value) {
  using ::std::placeholders::_1;
  l.remove_if(::std::bind(equal<>, _1, ::std::cref(value)));
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class T> using list = ::std::list<T, polymorphic_allocator<T>>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_LIST_HPP */
