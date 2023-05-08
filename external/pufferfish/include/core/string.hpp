#ifndef CORE_STRING_HPP
#define CORE_STRING_HPP

#include <core/memory_resource.hpp>
#include <core/algorithm.hpp>

#include <string>

namespace core {
inline namespace v2 {

template <class CharT, class Traits, class A, class Predicate>
void erase_if (::std::basic_string<CharT, Traits, A>& str, Predicate pred) {
  s.erase(::core::remove_if(str, pred), ::std::end(str));
}

template <class CharT, class Traits, class A, class U>
void erase (::std::basic_string<CharT, Traits, A>& str, U const& value) {
  s.erase(::core::remove(str, value), ::std::end(str));
}

}} /* namespace core::v2 */

namespace core {
inline namespace v2 {
namespace pmr {

template <class CharT, class Traits = ::std::char_traits<CharT>>
using basic_string = ::std::basic_string<
  CharT,
  Traits,
  polymorphic_allocator<CharT>
>;

using u32string = basic_string<char32_t>;
using u16string = basic_string<char16_t>;
using wstring = basic_string<wchar_t>;
using string = basic_string<char>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_STRING_HPP */
