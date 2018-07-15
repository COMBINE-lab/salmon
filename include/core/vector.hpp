#ifndef CORE_VECTOR_HPP
#define CORE_VECTOR_HPP

#include <core/memory_resource.hpp>
#include <core/algorithm.hpp>

#include <vector>

namespace core {
inline namespace v2 {

template <class T, class A, class Predicate>
void erase_if (::std::vector<T, A>& v, Predicate pred) {
  v.erase(::core::remove_if(v, pred), end(v));
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class T> using vector = ::std::vector<T, polymorphic_allocator<T>>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_VECTOR_HPP */
