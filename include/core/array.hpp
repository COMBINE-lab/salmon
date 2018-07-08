#ifndef CORE_ARRAY_HPP
#define CORE_ARRAY_HPP

#include <array>

#include <core/type_traits.hpp>
#include <core/functional.hpp>
#include <core/utility.hpp>
#include <core/meta.hpp>

namespace core {
inline namespace v2 {

template <class V = void, class... Args>
constexpr auto make_array (Args&&... args) -> ::std::array<
  meta::either<
    meta::all<
      ::std::is_void<V>::value,
      meta::none_of<meta::list<Args...>, is_reference_wrapper>()
    >(),
    common_type_t<Args...>,
    V
  >,
  sizeof...(Args)
> { return {{ core::forward<Args>(args)... }}; }

template <class T, ::std::size_t N, ::std::size_t... Is>
constexpr auto to_array (T (&array)[N], index_sequence<Is...>) -> ::std::array<
  remove_cv_t<T>, N
> { return {{ array[Is]... }}; }

template <class T, ::std::size_t N>
constexpr auto to_array (T (&array)[N]) -> ::std::array<remove_cv_t<T>, N> {
  return (to_array)(array, make_index_sequence<N> { });
}

}} /* namespace core::v2 */

#endif /* CORE_ARRAY_HPP */
