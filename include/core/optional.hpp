#ifndef CORE_OPTIONAL_HPP
#define CORE_OPTIONAL_HPP

#include <initializer_list>
#include <system_error>
#include <functional>
#include <memory>

#include <cstdlib>
#include <cstdint>

#include <core/type_traits.hpp>
#include <core/functional.hpp>
#include <core/utility.hpp>

#ifndef CORE_NO_EXCEPTIONS
#include <exception>
#include <stdexcept>
#endif /* CORE_NO_EXCEPTIONS */

#if defined(_MSC_VER)
  #pragma warning(push)
  #pragma warning(disable:4458)
  #pragma warning(disable:4702)
#endif /* defined(_MSC_VER) */

namespace core {
inline namespace v2 {
namespace impl {

template <class T>
using addressof_builtin_t = decltype(&::std::declval<T>());

template <class T>
using addressof_member_t = decltype(::std::declval<T>().operator&());

template <class T>
using addressof_free_t = decltype(operator &(::std::declval<T>()));

template <class T>
using has_no_addressof_overload = meta::all_t<
  meta::none<
    is_detected<addressof_member_t, T>::value,
    is_detected<addressof_free_t, T>::value
  >(),
  is_detected<addressof_builtin_t, T>::value
>;

struct place_t { };
constexpr place_t place { };

/* this is the default 'false' case */
template <class T, bool = ::std::is_trivially_destructible<T>::value>
struct storage {
  using value_type = T;
  static constexpr bool nothrow = ::std::is_nothrow_move_constructible<
    value_type
  >::value;

  union {
    ::std::uint8_t dummy;
    value_type val;
  };
  bool engaged { false };

  constexpr storage () noexcept : dummy { } { }
  storage (storage const& that) :
    engaged { that.engaged }
  {
    if (not this->engaged) { return; }
    ::new (::core::as_void(this->val)) value_type(that.val);
  }

  storage (storage&& that) noexcept(nothrow) :
    engaged { that.engaged }
  {
    if (not this->engaged) { return; }
    ::new (::core::as_void(this->val)) value_type(::core::move(that.val));
  }

  constexpr storage (value_type const& value) :
    val(value),
    engaged { true }
  { }

  constexpr storage (value_type&& value) noexcept(nothrow) :
    val(::core::move(value)),
    engaged { true }
  { }

  template <class... Args>
  constexpr explicit storage (place_t, Args&&... args) :
    val(::core::forward<Args>(args)...),
    engaged { true }
  { }

  ~storage () noexcept { if (this->engaged) { this->val.~value_type(); } }
};

template <class T>
struct storage<T, true> {
  using value_type = T;
  static constexpr bool nothrow = ::std::is_nothrow_move_constructible<
    value_type
  >::value;
  union {
    ::std::uint8_t dummy;
    value_type val;
  };
  bool engaged { false };

  constexpr storage () noexcept : dummy { } { }
  storage (storage const& that) :
    engaged { that.engaged }
  {
    if (not this->engaged) { return; }
    ::new (::core::as_void(this->val)) value_type(that.val);
  }

  storage (storage&& that) noexcept(nothrow) :
    engaged { that.engaged }
  {
    if (not this->engaged) { return; }
    ::new (::core::as_void(this->val)) value_type(::core::move(that.val));
  }

  constexpr storage (value_type const& value) :
    val(value),
    engaged { true }
  { }

  constexpr storage (value_type&& value) noexcept(nothrow) :
    val(::core::move(value)),
    engaged { true }
  { }

  template <class... Args>
  constexpr explicit storage (place_t, Args&&... args) :
    val(::core::forward<Args>(args)...),
    engaged { true }
  { }
};

} /* namespace impl */

struct in_place_t { };
struct nullopt_t { constexpr explicit nullopt_t (int) noexcept { } };

constexpr in_place_t in_place { };
constexpr nullopt_t nullopt { 0 };

#ifndef CORE_NO_EXCEPTIONS
struct bad_optional_access final : ::std::logic_error {
  using ::std::logic_error::logic_error;
};

struct bad_expected_type : ::std::logic_error {
  using ::std::logic_error::logic_error;
};

struct bad_result_condition final : ::std::logic_error {
  using ::std::logic_error::logic_error;
};

[[noreturn]] inline void throw_bad_optional_access () {
  throw bad_optional_access { "optional is disengaged" };
}

[[noreturn]] inline void throw_bad_result_condition () {
  throw bad_result_condition { "result<T> is valid" };
}

[[noreturn]] inline void throw_bad_void_result_condition () {
  throw bad_result_condition { "result<void> is valid" };
}

[[noreturn]] inline void throw_system_error (::std::error_condition e) {
  throw ::std::system_error { e.value(), e.category() };
}
#else /* CORE_NO_EXCEPTIONS */
[[noreturn]] inline void throw_bad_optional_access () { ::std::abort(); }
[[noreturn]] inline void throw_bad_result_condition () { ::std::abort(); }
[[noreturn]] inline void throw_bad_void_result_condition () { ::std::abort(); }
[[noreturn]] inline void throw_system_error (::std::error_condition) {
  ::std::abort();
}
#endif /* CORE_NO_EXCEPTIONS */

template <class Type>
struct optional final : private impl::storage<Type> {
  using base = impl::storage<Type>;
  using value_type = typename impl::storage<Type>::value_type;

  /* compiler enforcement */
  static_assert(
    not ::std::is_reference<value_type>::value,
    "Cannot have optional reference (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, nullopt_t>::value,
    "Cannot have optional<nullopt_t> (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, in_place_t>::value,
    "Cannot have optional<in_place_t> (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, ::std::nullptr_t>::value,
    "Cannot have optional nullptr (tautological)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, void>::value,
    "Cannot have optional<void> (ill-formed)"
  );

  static_assert(
    ::std::is_object<value_type>::value,
    "Cannot have optional with a non-object type (undefined behavior)"
  );

  static_assert(
    ::std::is_nothrow_destructible<value_type>::value,
    "Cannot have optional with non-noexcept destructible (undefined behavior)"
  );

  constexpr optional () noexcept { }
  optional (optional const&) = default;
  optional (optional&&) = default;
  ~optional () = default;

  constexpr optional (nullopt_t) noexcept { }

  constexpr optional (value_type const& value) :
    base { value }
  { }

  constexpr optional (value_type&& value) noexcept(base::nothrow) :
    base { ::core::move(value) }
  { }

  template <
    class... Args,
    class=enable_if_t< ::std::is_constructible<value_type, Args...>::value>
  > constexpr explicit optional (in_place_t, Args&&... args) :
    base { impl::place, ::core::forward<Args>(args)... }
  { }

  template <
    class T,
    class... Args,
    class=enable_if_t<
      ::std::is_constructible<
        value_type,
        ::std::initializer_list<T>&,
        Args...
      >::value
    >
  > constexpr explicit optional (
    in_place_t,
    ::std::initializer_list<T> il,
    Args&&... args
  ) : base { impl::place, il, ::core::forward<Args>(args)... } { }

  optional& operator = (optional const& that) {
    optional { that }.swap(*this);
    return *this;
  }

  optional& operator = (optional&& that) noexcept (
    meta::all<
      ::std::is_nothrow_move_assignable<value_type>::value,
      ::std::is_nothrow_move_constructible<value_type>::value
    >()
  ) {
    optional { ::core::move(that) }.swap(*this);
    return *this;
  }

  template <
    class T,
    class=enable_if_t<
      not ::std::is_same<decay_t<T>, optional>::value and
      ::std::is_constructible<value_type, T>::value and
      ::std::is_assignable<value_type&, T>::value
    >
  > optional& operator = (T&& value) {
    if (not this->engaged) { this->emplace(::core::forward<T>(value)); }
    else { **this = ::core::forward<T>(value); }
    return *this;
  }

  optional& operator = (nullopt_t) noexcept {
    if (this->engaged) {
      this->val.~value_type();
      this->engaged = false;
    }
    return *this;
  }

  void swap (optional& that) noexcept(
    meta::all<
      is_nothrow_swappable<value_type>::value,
      ::std::is_nothrow_move_constructible<value_type>::value
    >()
  ) {
    using ::std::swap;
    if (not *this and not that) { return; }
    if (*this and that) {
      swap(**this, *that);
      return;
    }

    auto& to_disengage = *this ? *this : that;
    auto& to_engage = *this ? that : *this;

    to_engage.emplace(::core::move(*to_disengage));
    to_disengage = nullopt;
  }

  constexpr explicit operator bool () const { return this->engaged; }

  constexpr value_type const& operator * () const noexcept {
    return this->val;
  }

  value_type& operator * () noexcept { return this->val; }

  constexpr value_type const* operator -> () const noexcept {
    return this->ptr(impl::has_no_addressof_overload<value_type> { });
  }

  value_type* operator -> () noexcept { return ::std::addressof(this->val); }

  template <class T, class... Args>
  void emplace (::std::initializer_list<T> il, Args&&... args) {
    *this = nullopt;
    ::new (::core::as_void(this->val)) value_type(
      il,
      ::core::forward<Args>(args)...
    );
    this->engaged = true;
  }

  template <class... Args>
  void emplace (Args&&... args) {
    *this = nullopt;
    ::new (::core::as_void(this->val)) value_type(
      ::core::forward<Args>(args)...
    );
    this->engaged = true;
  }

  constexpr value_type const& value () const noexcept(false) {
    return *this
      ? **this
      : (throw_bad_optional_access(), **this);
  }

  value_type& value () noexcept(false) {
    if (*this) { return **this; }
    throw_bad_optional_access();
  }

  template <
    class T,
    class=meta::when<
      meta::all<
        ::std::is_copy_constructible<value_type>::value,
        ::std::is_convertible<T, value_type>::value
      >()
    >
  > constexpr value_type value_or (T&& val) const& {
    return *this ? **this : static_cast<value_type>(::core::forward<T>(val));
  }

  template <
    class T,
    class=enable_if_t<
      meta::all<
        ::std::is_move_constructible<value_type>,
        ::std::is_convertible<T, value_type>
      >()
    >
  > value_type value_or (T&& val) && {
    return *this
      ? value_type { ::core::move(**this) }
      : static_cast<value_type>(::core::forward<T>(val));
  }

private:
  constexpr value_type const* ptr (::std::false_type) const noexcept {
    return &this->val;
  }

  value_type const* ptr (::std::true_type) const noexcept {
    return ::std::addressof(this->val);
  }
};

#ifndef CORE_NO_EXCEPTIONS
template <class Type>
struct expected final {
  using value_type = Type;
  using error_type = ::std::exception_ptr;

  static constexpr bool nothrow = ::std::is_nothrow_move_constructible<
    value_type
  >::value;

  /* compiler enforcement */
  static_assert(
    not ::std::is_reference<value_type>::value,
    "Cannot have expected reference (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, nullopt_t>::value,
    "Cannot have expected<nullopt_t> (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, in_place_t>::value,
    "Cannot have expected<in_place_t> (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, ::std::exception_ptr>::value,
    "Cannot have expected<std::exception_ptr> (tautological)"
  );

  static_assert(
    ::std::is_object<value_type>::value,
    "Cannot have expected with non-object type (undefined behavior)"
  );

  static_assert(
    ::std::is_nothrow_destructible<value_type>::value,
    "Cannot have expected with throwable destructor (undefined behavior)"
  );

  expected (::std::exception_ptr ptr) noexcept :
    ptr { ptr }
  { }

  expected (value_type const& val) :
    val(val),
    valid { true }
  { }

  expected (value_type&& val) noexcept(nothrow) :
    val(::core::move(val)),
    valid { true }
  { }

  template <
    class... Args,
    class=enable_if_t<::std::is_constructible<value_type, Args...>::value>
  > explicit expected (in_place_t, Args&&... args) :
    val(::core::forward<Args>(args)...),
    valid { true }
  { }

  template <
    class T,
    class... Args,
    class=enable_if_t<
      ::std::is_constructible<
        value_type,
        ::std::initializer_list<T>&,
        Args...
      >::value
    >
  > explicit expected (
    in_place_t,
    ::std::initializer_list<T> il,
    Args&&... args
  ) : val(il, ::core::forward<Args>(args)...), valid { true } { }

  expected (expected const& that) :
    valid { that.valid }
  {
    if (*this) { ::new (::core::as_void(this->val)) value_type(that.val); }
    else { ::new (::core::as_void(this->ptr)) error_type(that.ptr); }
  }

  expected (expected&& that) noexcept(nothrow) :
    valid { that.valid }
  {
    if (*this) {
      ::new (::core::as_void(this->val)) value_type(::core::move(that.val));
    } else { ::new (::core::as_void(this->ptr)) error_type(that.ptr); }
  }

  expected () :
    val { },
    valid { true }
  { }

  ~expected () noexcept { this->reset(); }

  expected& operator = (expected const& that) {
    expected { that }.swap(*this);
    return *this;
  }

  expected& operator = (expected&& that) noexcept(
    meta::all<
      ::std::is_nothrow_move_assignable<value_type>::value,
      ::std::is_nothrow_move_constructible<value_type>::value
    >()
  ) {
    expected { ::core::move(that) }.swap(*this);
    return *this;
  }

  template <
    class T,
    class=enable_if_t<
      meta::all<
        meta::none<
          meta::list<expected, ::std::exception_ptr>,
          ::std::is_same,
          decay_t<T>
        >(),
        ::std::is_constructible<value_type, T>::value,
        ::std::is_assignable<value_type&, T>::value
      >()
    >
  > expected& operator = (T&& value) {
    if (not *this) { this->emplace(::core::forward<T>(value)); }
    else { **this = ::core::forward<T>(value); }
    return *this;
  }

  expected& operator = (value_type const& value) {
    if (not *this) { this->emplace(value); }
    else { **this = value; }
    return *this;
  }

  expected& operator = (value_type&& value) {
    if (not *this) { this->emplace(::core::move(value)); }
    else { **this = ::core::move(value); }
    return *this;
  }

  expected& operator = (::std::exception_ptr ptr) {
    if (not *this) { this->ptr = ptr; }
    else {
      this->val.~value_type();
      ::new (::core::as_void(this->ptr)) ::std::exception_ptr(ptr);
      this->valid = false;
    }
    return *this;
  }

  void swap (expected& that) noexcept(
    meta::all<
      is_nothrow_swappable<value_type>::value,
      ::std::is_nothrow_move_constructible<value_type>::value
    >()
  ) {
    using ::std::swap;

    if (not *this and not that) {
      swap(this->ptr, that.ptr);
      return;
    }

    if (*this and that) {
      swap(this->val, that.val);
      return;
    }

    auto& to_invalidate = *this ? *this : that;
    auto& to_validate = *this ? that : *this;
    auto ptr = to_validate.ptr;
    to_validate.emplace(::core::move(*to_invalidate));
    to_invalidate = ptr;
  }

  explicit operator bool () const noexcept { return this->valid; }

  value_type const& operator * () const noexcept { return this->val; }
  value_type& operator * () noexcept { return this->val; }

  value_type const* operator -> () const noexcept {
    return ::std::addressof(this->val);
  }
  value_type* operator -> () noexcept { return ::std::addressof(this->val); }

  template <class T, class... Args>
  void emplace (::std::initializer_list<T> il, Args&&... args) {
    this->reset();
    ::new (::core::as_void(this->val)) value_type(
      il,
      ::core::forward<Args>(args)...
    );
    this->valid = true;
  }

  template <class... Args>
  void emplace (Args&&... args) {
    this->reset();
    ::new (::core::as_void(this->val)) value_type(
      ::core::forward<Args>(args)...
    );
    this->valid = true;
  }

  value_type const& value () const noexcept(false) {
    if (not *this) { ::std::rethrow_exception(this->ptr); }
    return **this;
  }

  value_type& value () noexcept(false) {
    if (not *this) { ::std::rethrow_exception(this->ptr); }
    return **this;
  }

  template <
    class T,
    class=meta::when<
      meta::all<
        ::std::is_copy_constructible<value_type>::value,
        ::std::is_convertible<T, value_type>::value
      >()
    >
  > value_type value_or (T&& val) const& {
    return *this ? **this : static_cast<value_type>(::core::forward<T>(val));
  }

  template <
    class T,
    class=enable_if_t<
      meta::all<
        ::std::is_move_constructible<value_type>::value,
        ::std::is_convertible<T, value_type>::value
      >()
    >
  > value_type value_or (T&& val) && {
    return *this
      ? value_type { ::core::move(**this) }
      : static_cast<value_type>(::core::forward<T>(val));
  }

  template <class E>
  E expect () const noexcept(false) {
    try { this->raise(); }
    catch (E const& e) { return e; }
    catch (...) {
      ::std::throw_with_nested(bad_expected_type { "unexpected exception" });
    }
  }

  [[noreturn]] void raise () const noexcept(false) {
    if (*this) { throw bad_expected_type { "expected<T> is valid" }; }
    ::std::rethrow_exception(this->ptr);
  }

  ::std::exception_ptr pointer () const noexcept(false) {
    if (*this) { throw bad_expected_type { "expected<T> is valid" }; }
    return this->ptr;
  }

private:

  void reset () {
    *this ? this->val.~value_type() : this->ptr.~exception_ptr();
  }

  union {
    value_type val;
    ::std::exception_ptr ptr;
  };
  bool valid { false };
};
#endif /* CORE_NO_EXCEPTIONS */


template <class Type>
struct result final {
  template <class U> friend struct result;

  using value_type = Type;
  using error_type = ::std::error_condition;

  static constexpr bool nothrow = ::std::is_nothrow_move_constructible<
    value_type
  >::value;

  static_assert(
    not ::std::is_same<decay_t<value_type>, nullopt_t>::value,
    "Cannot have result<nullopt_t> (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, in_place_t>::value,
    "Cannot have result<in_place_t> (ill-formed)"
  );

  static_assert(
    not ::std::is_same<decay_t<value_type>, ::std::error_condition>::value,
    "Cannot have result<error_condition> (tautological)"
  );

  static_assert(
    not ::std::is_error_condition_enum<value_type>::value,
    "Cannot have result with an error condition enum as the value (ill-formed)"
  );

  static_assert(
    not ::std::is_reference<value_type>::value,
    "Cannot have result<T&> (ill-formed)"
  );

  static_assert(
    ::std::is_object<value_type>::value,
    "Cannot have result with non-object type (undefined behavior)"
  );

  static_assert(
    ::std::is_nothrow_destructible<value_type>::value,
    "Cannot have result with throwable destructor (undefined behavior)"
  );

  result (int val, ::std::error_category const& cat) noexcept :
    valid { val == 0 }
  {
    if (*this) { ::new (::core::as_void(this->val)) value_type(); }
    else { ::new (::core::as_void(this->cnd)) error_type(val, cat); }
  }

  template <
    class E,
    class=meta::when<::std::is_error_condition_enum<E>::value>
  > result (E e) noexcept :
    valid { core::as_under(e) == 0 }
  {
    if (*this) { ::new (as_void(this->val)) value_type(); }
    else { ::new (as_void(this->cnd)) error_type(e); }
  }

  result (::std::error_condition const& ec) :
    valid { not ec }
  {
    if (*this) { ::new (::core::as_void(this->val)) value_type(); }
    else { ::new (::core::as_void(this->cnd)) error_type(ec); }
  }

  result (value_type const& val) :
    val(val),
    valid { true }
  { }

  result (value_type&& val) noexcept(nothrow) :
    val(::core::move(val)),
    valid { true }
  { }

  template <
    class T,
    class=meta::when<
      meta::all<
        meta::none<::std::is_same<result, result<T>>::value>(),
        ::std::is_constructible<
          value_type,
          add_lvalue_reference_t<add_const_t<T>>
        >::value
      >()
    >
  > result (result<T> const& that) :
    valid { that.valid }
  {
    if (*this) { ::new (::core::as_void(this->val)) value_type(that.val); }
    else { ::new (::core::as_void(this->cnd)) error_type(that.cnd); }
  }


  template <
    class T,
    class=meta::when<
      meta::all<
        meta::none<::std::is_same<result, result<T>>::value>(),
        ::std::is_constructible<value_type, add_rvalue_reference_t<T>>::value
      >()
    >
  > result (result<T>&& that) :
    valid { that.valid }
  {
    if (*this) {
      ::new (::core::as_void(this->val)) value_type(::core::move(that.val));
    } else { ::new (::core::as_void(this->cnd)) error_type(that.cnd); }
  }

  template <
    class... Args,
    class=enable_if_t<::std::is_constructible<value_type, Args...>::value>
  > explicit result (in_place_t, Args&&... args) :
    val(::core::forward<Args>(args)...),
    valid { true }
  { }

  template <
    class T,
    class... Args,
    class=enable_if_t<
      ::std::is_constructible<
        value_type,
        ::std::initializer_list<T>&,
        Args...
      >::value
    >
  > explicit result (
    in_place_t,
    ::std::initializer_list<T> il,
    Args&&... args
  ) : val(il, ::core::forward<Args>(args)...), valid { true } { }

  result (result const& that) :
    valid { that.valid }
  {
    if (*this) { ::new (::core::as_void(this->val)) value_type(that.val); }
    else { ::new (::core::as_void(this->cnd)) error_type(that.cnd); }
  }

  result (result&& that) :
    valid { that.valid }
  {
    if (*this) {
      ::new (::core::as_void(this->val)) value_type(::core::move(that.val));
    } else { ::new (::core::as_void(this->cnd)) error_type(that.cnd); }
  }

  result () :
    val { },
    valid { true }
  { }

  ~result () noexcept { this->reset(); }

  result& operator = (result const& that) {
    result { that }.swap(*this);
    return *this;
  }

  result& operator = (result&& that) noexcept(
    meta::all<
      ::std::is_nothrow_move_assignable<value_type>::value,
      ::std::is_nothrow_move_constructible<value_type>::value
    >()
  ) {
    result { ::core::move(that) }.swap(*this);
    return *this;
  }

  template <
    class T,
    class=meta::when<
      meta::all<
        meta::none_of<
          meta::list<result, ::std::error_condition>,
          ::std::is_same,
          decay_t<T>
        >(),
        ::std::is_constructible<value_type, T>::value,
        ::std::is_assignable<value_type&, T>::value
      >()
    >
  > result& operator = (T&& value) {
    if (not *this) { this->emplace(::core::forward<T>(value)); }
    else { **this = ::core::forward<T>(value); }
    return *this;
  }

  template <
    class ErrorConditionEnum,
    class=meta::when<::std::is_error_condition_enum<ErrorConditionEnum>::value>
  > result& operator = (ErrorConditionEnum e) {
    result { e }.swap(*this);
    return *this;
  }

  result& operator = (value_type const& value) {
    if (not *this) { this->emplace(value); }
    else { **this = value; }
    return *this;
  }

  result& operator = (value_type&& value) {
    if (not *this) { this->emplace(::core::move(value)); }
    else { **this = ::core::move(value); }
    return *this;
  }

  result& operator = (::std::error_condition const& cnd) {
    if (not cnd) { return *this; }
    if (not *this) { this->cnd = cnd; }
    else {
      this->reset();
      ::new (::core::as_void(this->cnd)) error_type(cnd);
      this->valid = false;
    }
    return *this;
  }

  void swap (result& that) noexcept(
    meta::all<
      is_nothrow_swappable<value_type>::value,
      ::std::is_nothrow_move_constructible<value_type>::value
    >()
  ) {
    using ::std::swap;
    if (not *this and not that) {
      swap(this->cnd, that.cnd);
      return;
    }

    if (*this and that) {
      swap(this->val, that.val);
      return;
    }

    auto& to_invalidate = *this ? *this : that;
    auto& to_validate = *this ? that : *this;
    auto cnd = to_validate.cnd;
    to_validate.emplace(::core::move(*to_invalidate));
    to_invalidate = cnd;
  }

  explicit operator bool () const noexcept { return this->valid; }

  value_type const& operator * () const noexcept { return this->val; }
  value_type& operator * () noexcept { return this->val; }

  value_type const* operator -> () const noexcept {
    return ::std::addressof(this->val);
  }

  value_type* operator -> () noexcept { return ::std::addressof(this->val); }

  template <class T, class... Args>
  void emplace (::std::initializer_list<T> il, Args&&... args) {
    this->reset();
    ::new (::core::as_void(this->val)) value_type(
      il,
      ::core::forward<Args>(args)...
    );
    this->valid = true;
  }

  template <class... Args>
  void emplace (Args&&... args) {
    this->reset();
    ::new (::core::as_void(this->val)) value_type(
      ::core::forward<Args>(args)...
    );
    this->valid = true;
  }

  value_type const& value () const noexcept(false) {
    if (*this) { return **this; }
    throw_system_error(this->cnd);
  }

  value_type& value () noexcept(false) {
    if (*this) { return **this; }
   throw_system_error(this->cnd);
  }

  template <
    class T,
    class=meta::when<
      meta::all<
        ::std::is_copy_constructible<value_type>::value,
        ::std::is_convertible<T, value_type>::value
      >()
    >
  > value_type value_or (T&& val) const& {
    return *this ? **this : static_cast<value_type>(::core::forward<T>(val));
  }

  template <
    class T,
    class=meta::when<
      meta::all<
        ::std::is_move_constructible<value_type>::value,
        ::std::is_convertible<T, value_type>::value
      >()
    >
  > value_type value_or (T&& val) && {
    return *this
      ? value_type { ::core::move(**this) }
      : static_cast<value_type>(::core::forward<T>(val));
  }

  ::std::error_condition const& condition () const noexcept(false) {
    if (*this) { throw_bad_result_condition(); }
    return this->cnd;
  }

private:
  void reset () {
    *this ? this->val.~value_type() : this->cnd.~error_condition();
  }

  union {
    value_type val;
    error_type cnd;
  };
  bool valid { false };
};

#ifndef CORE_NO_EXCEPTIONS
template <>
struct expected<void> final {
  using value_type = void;

  explicit expected (::std::exception_ptr ptr) noexcept : ptr { ptr } { }
  expected (expected const&) = default;
  expected (expected&&) = default;
  expected () = default;
  ~expected () = default;

  expected& operator = (::std::exception_ptr ptr) noexcept {
    expected { ptr }.swap(*this);
    return *this;
  }

  expected& operator = (expected const&) = default;
  expected& operator = (expected&&) = default;

  void swap (expected& that) noexcept {
    using ::std::swap;
    swap(this->ptr, that.ptr);
  }

  explicit operator bool () const noexcept { return not this->ptr; }

  template <class E>
  E expect () const noexcept(false) {
    try { this->raise(); }
    catch (E const& e) { return e; }
    catch (...) {
      ::std::throw_with_nested(bad_expected_type { "unexpected exception" });
    }
  }

  [[noreturn]] void raise () const noexcept(false) {
    if (*this) { throw bad_expected_type { "valid expected<void>" }; }
    ::std::rethrow_exception(this->ptr);
  }

  ::std::exception_ptr pointer () const noexcept(false) {
    if (*this) { throw bad_expected_type { "valid expected<void>" }; }
    return this->ptr;
  }

private:
  ::std::exception_ptr ptr;
};
#endif /* CORE_NO_EXCEPTIONS */

template <>
struct result<void> final {
  using value_type = void;
  using error_type = ::std::error_condition;

  result (int val, ::std::error_category const& cat) :
    cnd { val, cat }
  { }

  result (::std::error_condition const& ec) :
    cnd { ec }
  { }

  template <
    class ErrorConditionEnum,
    class=meta::when<::std::is_error_condition_enum<ErrorConditionEnum>::value>
  > result (ErrorConditionEnum e) noexcept :
    cnd { e }
  { }

  result (result const&) = default;
  result (result&&) = default;
  result () = default;

  template <
    class ErrorConditionEnum,
    class=enable_if_t<std::is_error_condition_enum<ErrorConditionEnum>::value>
  > result& operator = (ErrorConditionEnum e) noexcept {
    result { e }.swap(*this);
    return *this;
  }

  result& operator = (::std::error_condition const& ec) {
    result { ec }.swap(*this);
    return *this;
  }

  result& operator = (result const&) = default;
  result& operator = (result&&) = default;

  void swap (result& that) noexcept {
    using ::std::swap;
    swap(this->cnd, that.cnd);
  }

  explicit operator bool () const noexcept { return not this->cnd; }

  ::std::error_condition const& condition () const noexcept(false) {
    if (*this) { throw_bad_void_result_condition(); }
    return this->cnd;
  }

private:
  ::std::error_condition cnd;
};

/* operator == */
template <class T>
constexpr bool operator == (
  optional<T> const& lhs,
  optional<T> const& rhs
) noexcept {
  return static_cast<bool>(lhs) != static_cast<bool>(rhs)
    ? false
    : (not lhs and not rhs) or *lhs == *rhs;
}

template <class T>
constexpr bool operator == (optional<T> const& lhs, nullopt_t) noexcept {
  return not lhs;
}

template <class T>
constexpr bool operator == (nullopt_t, optional<T> const& rhs) noexcept {
  return not rhs;
}

template <class T>
constexpr bool operator == (optional<T> const& opt, T const& value) noexcept {
  return opt and *opt == value;
}

template <class T>
constexpr bool operator == (T const& value, optional<T> const& opt) noexcept {
  return opt and value == *opt;
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
bool operator == (expected<T> const& lhs, expected<T> const& rhs) noexcept {
  if (lhs and rhs) { return *lhs == *rhs; }
  return not lhs and not rhs;
}

template <class T>
bool operator == (expected<T> const& lhs, ::std::exception_ptr) noexcept {
  return not lhs;
}

template <class T>
bool operator == (::std::exception_ptr, expected<T> const& rhs) noexcept {
  return not rhs;
}

template <class T>
bool operator == (expected<T> const& lhs, T const& rhs) noexcept {
  return lhs and *lhs == rhs;
}

template <class T>
bool operator == (T const& lhs, expected<T> const& rhs) noexcept {
  return rhs and lhs == *rhs;
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
bool operator == (result<T> const& lhs, result<T> const& rhs) noexcept {
  if (lhs and rhs) { return *lhs == *rhs; }
  if (not lhs and not rhs) { return lhs.condition() == rhs.condition(); }
  return false;
}

template <class T>
bool operator == (
  result<T> const& lhs,
  ::std::error_condition const& rhs
) noexcept { return not lhs and rhs and lhs.condition() == rhs; }

template <class T>
bool operator == (
  ::std::error_condition const& lhs,
  result<T> const& rhs
) noexcept { return lhs and not rhs and lhs == rhs.condition(); }

template <class T>
bool operator == (
  result<T> const& lhs,
  ::std::error_code const& rhs
) noexcept { return not lhs and rhs and lhs.condition() == rhs; }

template <class T>
bool operator == (
  ::std::error_code const& lhs,
  result<T> const& rhs
) noexcept { return lhs and not rhs and lhs == rhs.condition(); }

template <class T>
bool operator == (result<T> const& res, T const& value) noexcept {
  return res and *res == value;
}

template <class T>
bool operator == (T const& value, result<T> const& res) noexcept {
  return res and value == *res;
}

#ifndef CORE_NO_EXCEPTIONS
/* void specializations */
template <>
inline bool operator == <void> (
  expected<void> const& lhs,
  expected<void> const& rhs
) noexcept {
  if (not lhs and not rhs) { return lhs.pointer() == rhs.pointer(); }
  return lhs and rhs;
}
#endif /* CORE_NO_EXCEPTIONS */

template <>
inline bool operator == <void> (
  result<void> const& lhs,
  result<void> const& rhs
) noexcept {
  if (not lhs and not rhs) { return lhs.condition() == rhs.condition(); }
  return lhs and rhs;
}

/* operator < */
template <class T>
constexpr bool operator < (
  optional<T> const& lhs,
  optional<T> const& rhs
) noexcept {
  return static_cast<bool>(rhs) == false ? false : not lhs or *lhs < *rhs;
}

template <class T>
constexpr bool operator < (optional<T> const& lhs, nullopt_t) noexcept {
  return not lhs;
}

template <class T>
constexpr bool operator < (nullopt_t, optional<T> const& rhs) noexcept {
  return static_cast<bool>(rhs);
}

template <class T>
constexpr bool operator < (optional<T> const& opt, T const& value) noexcept {
  return not opt or *opt < value;
}

template <class T>
constexpr bool operator < (T const& value, optional<T> const& opt) noexcept {
  return opt and value < *opt;
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
bool operator < (expected<T> const& lhs, expected<T> const& rhs) noexcept {
  if (not rhs) { return false; }
  return not lhs or *lhs < *rhs;
}

template <class T>
bool operator < (expected<T> const& lhs, ::std::exception_ptr) noexcept {
  return not lhs;
}

template <class T>
bool operator < (::std::exception_ptr, expected<T> const& rhs) noexcept {
  return static_cast<bool>(rhs);
}

template <class T>
bool operator < (expected<T> const& exp, T const& value) noexcept {
  return not exp or *exp < value;
}

template <class T>
bool operator < (T const& value, expected<T> const& exp) noexcept {
  return exp and value < *exp;
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
bool operator < (result<T> const& lhs, result<T> const& rhs) noexcept {
  if (not rhs and not lhs) { return lhs.condition() < rhs.condition(); }
  if (lhs and rhs) { return *lhs < *rhs; }
  return static_cast<bool>(rhs);
}

template <class T>
bool operator < (
  result<T> const& lhs,
  ::std::error_condition const& rhs
) noexcept { return not lhs and lhs.condition() < rhs; }

template <class T>
bool operator < (
  ::std::error_condition const& lhs,
  result<T> const& rhs
) noexcept { return rhs or lhs < rhs.condition(); }

template <>
inline bool operator < <void> (
  result<void> const& lhs,
  result<void> const& rhs
) noexcept {
  if (not lhs and not rhs) { return lhs.condition() < rhs.condition(); }
  return static_cast<bool>(rhs);
}

template <class T>
bool operator < (result<T> const& res, T const& value) noexcept {
  return not res or *res < value;
}

template <class T>
bool operator < (T const& value, result<T> const& res) noexcept {
  return res and value < *res;
}

/* operator != */
template <class T>
constexpr bool operator != (
  optional<T> const& lhs,
  optional<T> const& rhs
) noexcept { return not (lhs == rhs); }

template <class T>
constexpr bool operator != (optional<T> const& lhs, nullopt_t) noexcept {
  return not (lhs == nullopt);
}

template <class T>
constexpr bool operator != (nullopt_t, optional<T> const& rhs) noexcept {
  return not (nullopt == rhs);
}

template <class T>
constexpr bool operator != (optional<T> const& opt, T const& value) noexcept {
  return not (opt == value);
}

template <class T>
constexpr bool operator != (T const& value, optional<T> const& opt) noexcept {
  return not (value == opt);
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
bool operator != (expected<T> const& lhs, expected<T> const& rhs) noexcept {
  return not (lhs == rhs);
}

template <class T>
bool operator != (expected<T> const& lhs, ::std::exception_ptr rhs) noexcept {
  return not (lhs == rhs);
}

template <class T>
bool operator != (::std::exception_ptr lhs, expected<T> const& rhs) noexcept {
  return not (lhs == rhs);
}

template <class T>
bool operator != (expected<T> const& exp, T const& value) noexcept {
  return not (exp == value);
}

template <class T>
bool operator != (T const& value, expected<T> const& exp) noexcept {
  return not (value == exp);
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
bool operator != (result<T> const& lhs, result<T> const& rhs) noexcept {
  return not (lhs == rhs);
}

template <class T>
bool operator != (
  result<T> const& lhs,
  ::std::error_condition const& rhs
) noexcept { return not (lhs == rhs); }

template <class T>
bool operator != (
  ::std::error_condition const& lhs,
  result<T> const& rhs
) noexcept { return not (lhs == rhs); }

template <class T>
bool operator != (
  result<T> const& lhs,
  ::std::error_code const& rhs
) noexcept { return not (lhs == rhs); }

template <class T>
bool operator != (
  ::std::error_code const& lhs,
  result<T> const& rhs
) noexcept { return not (lhs == rhs); }

template <class T>
bool operator != (result<T> const& res, T const& value) noexcept {
  return not (res == value);
}

template <class T>
bool operator != (T const& value, result<T> const& res) noexcept {
  return not (value == res);
}

#ifndef CORE_NO_EXCEPTIONS
template <>
inline bool operator != <void> (
  expected<void> const& lhs,
  expected<void> const& rhs
) noexcept { return not (lhs == rhs); }
#endif /* CORE_NO_EXCEPTIONS */

template <>
inline bool operator != <void> (
  result<void> const& lhs,
  result<void> const& rhs
) noexcept { return not (lhs == rhs); }

/* optional<T> operator >= */
template <class T>
constexpr bool operator >= (
  optional<T> const& lhs,
  optional<T> const& rhs
) noexcept { return not (lhs < rhs); }

template <class T>
constexpr bool operator >= (optional<T> const&, nullopt_t) noexcept {
  return true;
}

template <class T>
constexpr bool operator >= (nullopt_t, optional<T> const& opt) noexcept {
  return opt < nullopt;
}

template <class T>
constexpr bool operator >= (optional<T> const& opt, T const& value) noexcept {
  return not (opt < value);
}

template <class T>
constexpr bool operator >= (T const& value, optional<T> const& opt) noexcept {
  return not (value < opt);
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
bool operator >= (expected<T> const& lhs, expected<T> const& rhs) noexcept {
  return not (lhs < rhs);
}

template <class T>
bool operator >= (expected<T> const&, ::std::exception_ptr) noexcept {
  return true;
}

template <class T>
bool operator >= (::std::exception_ptr ptr, expected<T> const& exp) noexcept {
  return exp < ptr;
}

template <class T>
bool operator >= (expected<T> const& exp, T const& value) noexcept {
  return not (exp < value);
}

template <class T>
bool operator >= (T const& value, expected<T> const& exp) noexcept {
  return not (value < exp);
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
bool operator >= (result<T> const& lhs, result<T> const& rhs) noexcept {
  return not (lhs < rhs);
}

template <class T>
bool operator >= (
  result<T> const& lhs,
  ::std::error_condition const& rhs
) noexcept { return not (lhs < rhs); }

template <class T>
bool operator >= (
  ::std::error_condition const& lhs,
  result<T> const& rhs
) noexcept { return not (lhs < rhs); }

template <class T>
bool operator >= (result<T> const& res, T const& value) noexcept {
  return not (res < value);
}

template <class T>
bool operator >= (T const& value, result<T> const& res) noexcept {
  return not (value < res);
}

template <>
inline bool operator >= <void> (
  result<void> const& lhs,
  result<void> const& rhs
) noexcept { return not (lhs < rhs); }

/* operator <= */
template <class T>
constexpr bool operator <= (
  optional<T> const& lhs,
  optional<T> const& rhs
) noexcept { return not (rhs < lhs); }

template <class T>
constexpr bool operator <= (optional<T> const& lhs, nullopt_t) noexcept {
  return not lhs;
}

template <class T>
constexpr bool operator <= (nullopt_t, optional<T> const&) noexcept {
  return true;
}

template <class T>
constexpr bool operator <= (optional<T> const& opt, T const& value) noexcept {
  return not (opt > value);
}

template <class T>
constexpr bool operator <= (T const& value, optional<T> const& opt) noexcept {
  return not (value > opt);
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
bool operator <= (expected<T> const& lhs, expected<T> const& rhs) noexcept {
  return not (rhs < lhs);
}

template <class T>
bool operator <= (expected<T> const& lhs, ::std::exception_ptr) noexcept {
  return not lhs;
}

template <class T>
bool operator <= (::std::exception_ptr, expected<T> const&) noexcept {
  return true;
}

template <class T>
bool operator <= (expected<T> const& exp, T const& value) noexcept {
  return not (value < exp);
}

template <class T>
bool operator <= (T const& value, expected<T> const& exp) noexcept {
  return not (exp < value);
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
bool operator <= (result<T> const& lhs, result<T> const& rhs) noexcept {
  return not (rhs < lhs);
}

template <class T>
bool operator <= (
  result<T> const& lhs,
  ::std::error_condition const& rhs
) noexcept { return not (rhs < lhs); }

template <class T>
bool operator <= (
  ::std::error_condition const& lhs,
  result<T> const& rhs
) noexcept { return not (rhs < lhs); }

template <class T>
bool operator <= (result<T> const& res, T const& value) noexcept {
  return not (value < res);
}

template <class T>
bool operator <= (T const& value, result<T> const& res) noexcept {
  return not (res < value);
}

/* operator > */
template <class T>
constexpr bool operator > (
  optional<T> const& lhs,
  optional<T> const& rhs
) noexcept { return rhs < lhs; }

template <class T>
constexpr bool operator > (optional<T> const& lhs, nullopt_t) noexcept {
  return static_cast<bool>(lhs);
}

template <class T>
constexpr bool operator > (nullopt_t, optional<T> const&) noexcept {
  return false;
}

template <class T>
constexpr bool operator > (optional<T> const& opt, T const& value) noexcept {
  return value < opt;
}

template <class T>
constexpr bool operator > (T const& value, optional<T> const& opt) noexcept {
  return opt < value;
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
bool operator > (expected<T> const& lhs, expected<T> const& rhs) noexcept {
  return rhs < lhs;
}

template <class T>
bool operator > (expected<T> const& exp, T const& value) noexcept {
  return value < exp;
}

template <class T>
bool operator > (T const& value, expected<T> const& exp) noexcept {
  return exp < value;
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
bool operator > (result<T> const& lhs, result<T> const& rhs) noexcept {
  return rhs < lhs;
}

template <class T>
bool operator > (
  result<T> const& lhs,
  ::std::error_condition const& rhs
) noexcept { return rhs < lhs; }

template <class T>
bool operator > (
  ::std::error_condition const& lhs,
  result<T> const& rhs
) noexcept { return rhs < lhs; }

template <class T>
bool operator > (result<T> const& res, T const& value) noexcept {
  return value < res;
}

template <class T>
bool operator > (T const& value, result<T> const& res) noexcept {
  return res < value;
}

/* make_ functions */
template <class Type>
auto make_optional (Type&& value) -> optional<decay_t<Type>> {
  return optional<decay_t<Type>> { ::core::forward<Type>(value) };
}

#ifndef CORE_NO_EXCEPTIONS
template <class T>
auto make_expected (::std::exception_ptr error) -> expected<T> {
  return expected<T> { error };
}

template <class T>
auto make_expected (T&& value) -> enable_if_t<
  not ::std::is_base_of< ::std::exception, decay_t<T>>::value,
  expected<decay_t<T>>
> { return expected<T> { ::core::forward<T>(value) }; }

template <class T, class U>
auto make_expected (U&& value) -> enable_if_t<
  ::std::is_base_of< ::std::exception, decay_t<U>>::value,
  expected<T>
> {
  return make_expected<T>(
    ::std::make_exception_ptr(::core::forward<U>(value))
  );
}
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
auto make_result (::std::error_condition const& e) -> result<T> {
  return result<T> { e };
}

template <
  class T,
  class ErrorConditionEnum,
  class=enable_if_t<
    ::std::is_error_condition_enum<ErrorConditionEnum>::value
  >
> auto make_result (ErrorConditionEnum e) -> result<T> {
  return result<T> { e };
}

template <class T> auto make_result (T&& value) -> result<decay_t<T>> {
  return result<T> { ::core::forward<T>(value) };
}

template <class T>
result<T> make_result (int val, ::std::error_category const& cat) {
  return make_result<T>(::std::error_condition { val, cat });
}

template <class T>
void swap (optional<T>& lhs, optional<T>& rhs) noexcept(
  noexcept(lhs.swap(rhs))
) { lhs.swap(rhs); }

#ifndef CORE_NO_EXCEPTIONS
template <class T>
void swap (expected<T>& lhs, expected<T>& rhs) noexcept(
  noexcept(lhs.swap(rhs))
) { lhs.swap(rhs); }
#endif /* CORE_NO_EXCEPTIONS */

template <class T>
void swap (result<T>& lhs, result<T>& rhs) noexcept (
  noexcept(lhs.swap(rhs))
) { lhs.swap(rhs); }

}} /* namespace core::v2 */

namespace std {

template <class Type>
struct hash<::core::v2::optional<Type>> {
  using result_type = typename hash<Type>::result_type;
  using argument_type = ::core::v2::optional<Type>;

  result_type operator () (argument_type const& value) const noexcept {
    return value ? hash<Type> { }(*value) : result_type { };
  }
};

#ifndef CORE_NO_EXCEPTIONS
template <class Type>
struct hash<::core::v2::expected<Type>> {
  using result_type = typename hash<Type>::result_type;
  using argument_type = ::core::v2::expected<Type>;

  result_type operator () (argument_type const& value) const noexcept {
    return value ? hash<Type> { }(*value) : result_type { };
  }
};
#endif /* CORE_NO_EXCEPTIONS */

template <class Type>
struct hash<::core::v2::result<Type>> {
  using result_type = typename hash<Type>::result_type;
  using argument_type = ::core::v2::result<Type>;

  result_type operator () (argument_type const& value) const noexcept {
    return value ? hash<Type> { }(*value) : result_type { };
  }
};

} /* namespace std */

#if defined(_MSC_VER)
  #pragma warning(pop)
#endif /* defined(_MSC_VER) */

#endif /* CORE_OPTIONAL_HPP */
