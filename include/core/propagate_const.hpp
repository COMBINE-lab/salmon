#ifndef CORE_PROPAGATE_CONST_HPP
#define CORE_PROPAGATE_CONST_HPP

#include <core/type_traits.hpp>
#include <core/functional.hpp>
#include <core/utility.hpp>
#include <core/meta.hpp>

#include <memory>

namespace core {
inline namespace v2 {

namespace impl {

template <class T>
constexpr auto pointer_to (T const& t) -> decltype(t.get()) { return t.get(); }
template <class T> constexpr T* pointer_to (T* t) { return t; }

} /* namespace impl */

template <class T>
struct propagate_const final {

  using element_type = remove_reference_t<decltype(*::std::declval<T&>())>;
  using underlying_type = T;

  /* Although these are public, they are undocumented and for internal use */
  template <class U>
  using is_self = ::std::is_same<decay_t<U>, propagate_const>;
  using is_nothrow_swap = is_nothrow_swappable<underlying_type>;

  propagate_const (propagate_const const&) = delete;
  constexpr propagate_const (propagate_const&&) = default;
  constexpr propagate_const () = default;

  template <
    class U,
    meta::require<::std::is_constructible<T, U&&>::value> = __LINE__,
    meta::inhibit<::std::is_convertible<U&&, T>::value> = __LINE__
  > explicit propagate_const (propagate_const<U>&& that) :
    pointer { ::core::move(that.pointer) }
  { }

  template <
    class U,
    meta::require<not ::std::is_convertible<U&&, T>::value> = __LINE__,
    meta::require<::std::is_constructible<T, U&&>::value> = __LINE__,
    meta::require<not is_self<U>::value> = __LINE__
  > explicit propagate_const (U&& value) :
    pointer { ::core::forward<U>(value) }
  { }

  template <
    class U,
    meta::require<::std::is_constructible<T, U&&>::value> = __LINE__,
    meta::require<::std::is_convertible<U&&, T>::value> = __LINE__
  > propagate_const (propagate_const<U>&& that) :
    pointer { move(that.pointer) }
  { }

  template <
    class U,
    meta::require<::std::is_constructible<T, U&&>::value> = __LINE__,
    meta::require<::std::is_convertible<U&&, T>::value> = __LINE__,
    meta::require<not is_self<U>::value> = __LINE__
  > propagate_const (U&& value) :
    pointer { ::core::forward<U>(value) }
  { }

  propagate_const& operator = (propagate_const const&) = delete;
  propagate_const& operator = (propagate_const&&) = default;

  template <class U>
  propagate_const& operator = (propagate_const<U>&& that) {
    this->pointer = ::core::move(that.pointer);
    return *this;
  }

  template <
    class U,
    meta::require<::std::is_convertible<U, T>::value> = __LINE__,
    meta::require<not is_self<U>::value> = __LINE__
  > propagate_const& operator = (U&& value) {
    this->pointer = ::core::forward<U>(value);
    return *this;
  }

  void swap (propagate_const& that) noexcept(is_nothrow_swap::value) {
    using ::std::swap;
    swap(this->pointer, that.pointer);
  }

  constexpr explicit operator bool () const noexcept {
    return bool(this->pointer);
  }

  constexpr operator element_type const* () const { return this->get(); }
  operator element_type* () { return this->get(); }

  constexpr element_type const* operator -> () const { return this->get(); }
  element_type* operator -> () { return this->get(); }

  constexpr element_type const& operator * () const { return *this->get(); }
  element_type& operator * () { return *this->get(); }

  constexpr element_type const* get () const {
    return impl::pointer_to(this->pointer);
  }

  element_type* get () { return impl::pointer_to(this->pointer); }

  friend constexpr underlying_type const& get_underlying (
    propagate_const const& pc
  ) noexcept { return pc.pointer; }

  friend constexpr underlying_type& get_underlying (
    propagate_const& pc
  ) noexcept { return pc.pointer; }

private:
  underlying_type pointer;
};

template <class T>
void swap (propagate_const<T>& lhs, propagate_const<T>& rhs) noexcept(
  noexcept(lhs.swap(rhs))
) { lhs.swap(rhs); }

template <class T>
constexpr bool operator == (
  propagate_const<T> const& l,
  propagate_const<T> const& r
) { return get_underlying(l) == get_underlying(r); }

template <class T>
constexpr bool operator != (
  propagate_const<T> const& l,
  propagate_const<T> const& r
) { return get_underlying(l) != get_underlying(r); }

template <class T>
constexpr bool operator >= (
  propagate_const<T> const& l,
  propagate_const<T> const& r
) { return get_underlying(l) >= get_underlying(r); }

template <class T>
constexpr bool operator <= (
  propagate_const<T> const& l,
  propagate_const<T> const& r
) { return get_underlying(l) <= get_underlying(r); }

template <class T>
constexpr bool operator > (
  propagate_const<T> const& l,
  propagate_const<T> const& r
) { return get_underlying(l) > get_underlying(r); }

template <class T>
constexpr bool operator < (
  propagate_const<T> const& l,
  propagate_const<T> const& r
) { return get_underlying(l) < get_underlying(r); }

template <class T, class U>
constexpr bool operator == (propagate_const<T> const& l, U const& r) {
  return get_underlying(l) == r;
}

template <class T, class U>
constexpr bool operator != (propagate_const<T> const& l, U const& r) {
  return get_underlying(l) != r;
}

template <class T, class U>
constexpr bool operator >= (propagate_const<T> const& l, U const& r) {
  return get_underlying(l) >= r;
}

template <class T, class U>
constexpr bool operator <= (propagate_const<T> const& l, U const& r) {
  return get_underlying(l) <= r;
}

template <class T, class U>
constexpr bool operator > (propagate_const<T> const& l, U const& r) {
  return get_underlying(l) > r;
}

template <class T, class U>
constexpr bool operator < (propagate_const<T> const& l, U const& r) {
  return get_underlying(l) < r;
}

template <class T>
constexpr bool operator == (
  propagate_const<T> const& l,
  ::std::nullptr_t
) { return not l; }

template <class T>
constexpr bool operator != (
  propagate_const<T> const& l,
  ::std::nullptr_t
) { return bool(l); }

template <class T>
constexpr bool operator == (
  ::std::nullptr_t,
  propagate_const<T> const& r
) { return not r; }

template <class T>
struct equal_to<propagate_const<T>> {
  constexpr bool operator () (
    propagate_const<T> const& l,
    propagate_const<T> const& r
  ) const { return equal_to<T> { }(get_underlying(l), get_underlying(r)); }
};

template <class T>
struct not_equal_to<propagate_const<T>> {
  constexpr bool operator () (
    propagate_const<T> const& l,
    propagate_const<T> const& r
  ) const { return not_equal_to<T> { }(get_underlying(l), get_underlying(r)); }
};

template <class T>
struct greater_equal<propagate_const<T>> {
  constexpr bool operator () (
    propagate_const<T> const& l,
    propagate_const<T> const& r
  ) const {
    return greater_equal<T> { }(get_underlying(l), get_underlying(r));
  }
};

template <class T>
struct less_equal<propagate_const<T>> {
  constexpr bool operator () (
    propagate_const<T> const& l,
    propagate_const<T> const& r
  ) const { return less_equal<T> { }(get_underlying(l), get_underlying(r)); }
};

template <class T>
struct greater<propagate_const<T>> {
  constexpr bool operator () (
    propagate_const<T> const& l,
    propagate_const<T> const& r
  ) const { return greater<T> { }(get_underlying(l), get_underlying(r)); }
};

template <class T>
struct less<propagate_const<T>> {
  constexpr bool operator () (
    propagate_const<T> const& l,
    propagate_const<T> const& r
  ) const { return less<T> { }(get_underlying(l), get_underlying(r)); }
};

}} /* namespace core::v2 */

namespace std {

template <class T>
struct hash<::core::v2::propagate_const<T>> {
  using argument_type = ::core::v2::propagate_const<T>;
  using result_type = size_t;

  result_type operator () (argument_type const& value) const {
    using underlying_type = typename argument_type::underlying_type;
    return hash<underlying_type> { }(get_underlying(value));
  }
};

template <class T>
struct equal_to<::core::v2::propagate_const<T>> {
  using result_type = bool;
  using first_argument_type = ::core::v2::propagate_const<T>;
  using second_argument_type = ::core::v2::propagate_const<T>;

  bool operator () (
    ::core::v2::propagate_const<T> const& lhs,
    ::core::v2::propagate_const<T> const& rhs
  ) const { return equal_to<T> { }(get_underlying(lhs), get_underlying(rhs)); }
};

template <class T>
struct not_equal_to<::core::v2::propagate_const<T>> {
  using result_type = bool;
  using first_argument_type = ::core::v2::propagate_const<T>;
  using second_argument_type = ::core::v2::propagate_const<T>;

  bool operator () (
    ::core::v2::propagate_const<T> const& lhs,
    ::core::v2::propagate_const<T> const& rhs
  ) const {
    return not_equal_to<T> { }(get_underlying(lhs), get_underlying(rhs));
  }
};

template <class T>
struct greater_equal<::core::v2::propagate_const<T>> {
  using result_type = bool;
  using first_argument_type = ::core::v2::propagate_const<T>;
  using second_argument_type = ::core::v2::propagate_const<T>;

  bool operator () (
    ::core::v2::propagate_const<T> const& lhs,
    ::core::v2::propagate_const<T> const& rhs
  ) const {
    return greater_equal<T> { }(get_underlying(lhs), get_underlying(rhs));
  }
};

template <class T>
struct less_equal<::core::v2::propagate_const<T>> {
  using result_type = bool;
  using first_argument_type = ::core::v2::propagate_const<T>;
  using second_argument_type = ::core::v2::propagate_const<T>;

  bool operator () (
    ::core::v2::propagate_const<T> const& lhs,
    ::core::v2::propagate_const<T> const& rhs
  ) const {
    return less_equal<T> { }(get_underlying(lhs), get_underlying(rhs));
  }
};

template <class T>
struct greater<::core::v2::propagate_const<T>> {
  using result_type = bool;
  using first_argument_type = ::core::v2::propagate_const<T>;
  using second_argument_type = ::core::v2::propagate_const<T>;

  bool operator () (
    ::core::v2::propagate_const<T> const& lhs,
    ::core::v2::propagate_const<T> const& rhs
  ) const { return greater<T> { }(get_underlying(lhs), get_underlying(rhs)); }
};

template <class T>
struct less<::core::v2::propagate_const<T>> {
  using result_type = bool;
  using first_argument_type = ::core::v2::propagate_const<T>;
  using second_argument_type = ::core::v2::propagate_const<T>;

  bool operator () (
    ::core::v2::propagate_const<T> const& lhs,
    ::core::v2::propagate_const<T> const& rhs
  ) const { return less<T> { }(get_underlying(lhs), get_underlying(rhs)); }
};

} /* namespace std */

#endif /* CORE_PROPAGATE_CONST_HPP */
