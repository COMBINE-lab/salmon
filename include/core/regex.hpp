#ifndef CORE_REGEX_HPP
#define CORE_REGEX_HPP

#include <core/memory_resource.hpp>
#include <core/string.hpp>

#include <regex>

namespace core {
inline namespace v2 {
namespace pmr {

template <class BidirIt>
using match_results = ::std::match_results<
  BidirIt,
  polymorphic_allocator<::std:::sub_match<BidirIt>>
>;

using wsmatch = match_results<wstring::const_iterator>;
using wcmatch = match_results<wchar_t const*>;
using smatch = match_results<string::const_iterator>;
using cmatch = match_results<char const*>;

}}} /* namespace core::v2::pmr */

#endif /* CORE_REGEX_HPP */
