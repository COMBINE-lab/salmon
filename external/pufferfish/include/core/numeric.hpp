#ifndef CORE_NUMERIC_HPP
#define CORE_NUMERIC_HPP

#include <numeric>

#include <core/range.hpp>

namespace core {
inline namespace v2 {

template <class Range, class T>
auto iota (Range&& rng, T&& value) -> enable_if_t<is_range<Range>::value> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_forward = decltype(range)::is_forward;
  static_assert(is_forward, "iota requires ForwardIterators");
  return ::std::iota(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<T>(value)
  );
}

template <class Range, class T>
auto accumulate (Range&& rng, T&& init) -> enable_if_t<
  is_range<Range>::value,
  decay_t<T>
> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "accumulate requires InputIterators");
  return ::std::accumulate(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<T>(init)
  );
}

template <class Range, class T, class BinaryOp>
auto accumulate (Range&& rng, T&& init, BinaryOp&& op) -> enable_if_t<
  is_range<Range>::value,
  decay_t<T>
> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "accumulate requires InputIterators");
  return ::std::accumulate(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<T>(init),
    ::std::forward<BinaryOp>(op)
  );
}

template <class Range, class InputIt, class T>
auto inner_product (Range&& rng, InputIt&& it, T&& value) -> enable_if_t<
  is_range<Range>::value,
  decay_t<T>
> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "inner_product requires InputIterators");
  return ::std::inner_product(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<InputIt>(it),
    ::std::forward<T>(value)
  );
}

template <class Range, class InputIt, class T, class BinaryOp, class BinaryOp2>
auto inner_product (
  Range&& rng,
  InputIt&& it,
  T&& value,
  BinaryOp&& op,
  BinaryOp2&& op2
) -> enable_if_t<is_range<Range>::value, decay_t<T>> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "inner_product requires InputIterators");
  return ::std::inner_product(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<InputIt>(it),
    ::std::forward<T>(value),
    ::std::forward<BinaryOp>(op),
    ::std::forward<BinaryOp2>(op2)
  );
}

template <class Range, class OutputIt>
auto adjacent_difference (Range&& rng, OutputIt&& it) -> enable_if_t<
  is_range<Range>::value,
  decay_t<OutputIt>
> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "adjacent_difference requires InputIterators");
  return ::std::adjacent_difference(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<OutputIt>(it)
  );
}

template <class Range, class OutputIt, class BinaryOp>
auto adjacent_difference (
  Range&& rng,
  OutputIt&& it,
  BinaryOp&& op
) -> enable_if_t<is_range<Range>::value, decay_t<OutputIt>> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "adjacent_difference requires InputIterators");
  return ::std::adjacent_difference(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<OutputIt>(it),
    ::std::forward<BinaryOp>(op)
  );
}

template <class Range, class OutputIt>
auto partial_sum (Range&& rng, OutputIt&& it) -> enable_if_t<
  is_range<Range>::value,
  decay_t<OutputIt>
> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "partial_sum requires InputIterators");
  return ::std::partial_sum(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<OutputIt>(it)
  );
}

template <class Range, class OutputIt, class BinaryOp>
auto partial_sum (Range&& rng, OutputIt&& it, BinaryOp&& op) -> enable_if_t<
  is_range<Range>::value,
  decay_t<OutputIt>
> {
  auto range = make_range(::std::forward<Range>(rng));
  constexpr auto is_input = decltype(range)::is_input;
  static_assert(is_input, "partial_sum requires InputIterators");
  return ::std::partial_sum(
    ::std::begin(range),
    ::std::end(range),
    ::std::forward<OutputIt>(it),
    ::std::forward<BinaryOp>(op)
  );
}

/* Because they need to be constexpr, these are EXTREMELY ineffecient */
template <class M, class N>
constexpr auto gcd (M m, N n) -> meta::when<
  meta::all_of<meta::list<M, N>, ::std::is_integral>(),
  common_type_t<M, N>
> { return m % n ? gcd(m, m % n) : n; }

template <class M, class N>
constexpr auto lcm (M m, N n) -> meta::when<
  meta::all_of<meta::list<M, N>, std::is_integral>(),
  common_type_t<M, N>
> { return m / gcd(m, n) * n; }

}} /* namespace core::v2 */

#endif /* CORE_NUMERIC_HPP */
