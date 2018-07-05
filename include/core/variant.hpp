#ifndef CORE_VARIANT_HPP
#define CORE_VARIANT_HPP

#include <core/type_traits.hpp>
#include <core/functional.hpp>
#include <core/typeinfo.hpp>
#include <core/utility.hpp>
#include <core/array.hpp>

#include <stdexcept>
#include <typeinfo>
#include <limits>

#include <cstdlib>
#include <cstdint>

namespace core {
inline namespace v2 {
namespace impl {

/* This is used to get around GCC's inability to expand lambdas in variadic
 * template functions. You make me so sad sometimes, GCC.
 */
template <class Visitor, class Type, class Data, class Result, class... Args>
auto gen () -> Result {
  return [](Visitor&& visitor, Data* data, Args&&... args) {
    return ::core::v2::invoke(
      ::core::forward<Visitor>(visitor),
      *static_cast<Type*>(data),
      ::core::forward<Args>(args)...
    );
  };
}

template <class, class, class...> struct result;
template <class V, class... Ts, class... Args>
struct result<V, meta::list<Ts...>, Args...> final : ::std::conditional<
  meta::all<
    meta::all_of<
      meta::list<result_of_t<V(Ts, Args...)>...>,
      ::std::is_same,
      result_of_t<V(Ts, Args...)>
    >()...
  >(),
  result_of_t<V(meta::get<meta::list<Ts...>, 0>, Args...)>,
  common_type_t<result_of_t<V(Ts, Args...)>...>
> { };

template <class V, class T, class... Args>
using result_t = typename result<V, T, Args...>::type;

/* Used to provide lambda based 'pattern matching' for variant and optional
 * types.
 *
 * Based off of Dave Abrahams C++11 'generic lambda' example (no longer
 * available on the internet)
 */

template <class... Lambdas> struct overload;
template <class Lambda> struct overload<Lambda> : Lambda {
  using call_type = Lambda;
  using call_type::operator ();
};

template <class Lambda, class... Lambdas>
struct overload<Lambda, Lambdas...> :
  private Lambda,
  private overload<Lambdas...>::call_type
{
  using base_type = typename overload<Lambdas...>::call_type;

  using lambda_type = Lambda;
  using call_type = overload;

  overload (Lambda&& lambda, Lambdas&&... lambdas) :
    lambda_type(pass<Lambda>(lambda)),
    base_type(pass<Lambdas>(lambdas)...)
  { }

  using lambda_type::operator ();
  using base_type::operator ();
};

template <class... Lambdas>
auto make_overload(Lambdas&&... lambdas) -> overload<Lambdas...> {
  return overload<Lambdas...> { pass<Lambdas>(lambdas)... };
}

}}} /* namespace core::v2::impl */

namespace core {
inline namespace v2 {

#ifndef CORE_NO_EXCEPTIONS
struct bad_variant_get final : ::std::logic_error {
  using ::std::logic_error::logic_error;
};
[[noreturn]] inline void throw_bad_variant_get () {
  throw bad_variant_get { "incorrect type" };
}
#else /* CORE_NO_EXCEPTIONS */
[[noreturn]] inline void throw_bad_variant_get () { ::std::abort(); }
#endif /* CORE_NO_EXCEPTIONS */

template <::std::size_t> struct emplace_index_t { };
template <class T> struct emplace_type_t { };

/* visitation semantics require that, given a callable type C, and variadic
 * arguments Args... that the return type of the visit will be SFINAE-ified
 * as common_type_t<invoke_of_t<C, Args>...> (this assumes a variadic
 * approach can be taken with common_type, index it cannot at this time. A
 * custom SFINAE-capable version has been written within the type traits
 * component.
 *
 * Obviously if a common type cannot be found, then the visitation function
 * cannot be generated.
 *
 * These same semantics are required for variant<Ts...>::match index simply
 * calls visit with a generate overload<Lambdas...> type.
 */
template <class... Ts>
class variant final {
  using typelist = meta::list<Ts...>;

  //static_assert(meta::all<(meta::count<typelist, Ts>() == 1)...>(), "");
  static_assert(meta::none_of<typelist, ::std::is_reference>(), "");
  static_assert(meta::none_of<typelist, ::std::is_void>(), "");

  using storage_type = aligned_union_t<0, Ts...>;

  template <::std::size_t N> using size = meta::integral<::std::size_t, N>;
  template <::std::size_t N> using element = meta::get<typelist, N>;

  struct copier final {
    using data_type = add_pointer_t<void>;
    data_type data;

    template <class T>
    void operator ()(T const& value) const { ::new (this->data) T(value); }
  };

  struct mover final {
    using data_type = add_pointer_t<void>;
    data_type data;

    template <class T>
    void operator () (T&& value) {
      ::new (this->data) decay_t<T>(::core::move(value));
    }
  };

  struct destroyer final {
    template <class T> void operator ()(T const& value) const { value.~T(); }
  };

  struct swapper final {
    using data_type = add_pointer_t<void>;
    data_type data;

    template <class T> void operator () (T const&) = delete;

    template <class T>
    void operator ()(T& value) noexcept(is_nothrow_swappable<T>::value) {
      using ::std::swap;
      swap(*static_cast<T*>(this->data), value);
    }
  };

  struct equality final {
    using data_type = add_pointer_t<add_const_t<void>>;
    data_type data;

    template <class T>
    bool operator ()(T const& value) {
      return equal_to<> { }(*static_cast<T const*>(this->data), value);
    }
  };

  struct less_than final {
    using data_type = add_pointer_t<add_const_t<void>>;
    data_type data;

    template <class T>
    bool operator ()(T const& value) noexcept {
      return less<> { }(*static_cast<T const*>(this->data), value);
    }
  };

  struct typeinfo final {
    template <class T>
    type_info const* operator ()(T&&) const noexcept {
      return ::std::addressof(type_of<decay_t<T>>());
    }
  };

  template <
    ::std::size_t N,
    class=enable_if_t<N < typelist::size()>,
    class T
  > explicit variant (size<N>&&, ::std::false_type&&, T&& value) :
    variant {
      size<N + 1> { },
      ::std::is_constructible<meta::get<meta::list<Ts...>, N + 1>, T> { },
      ::core::forward<T>(value)
    }
  { }

  template <
    ::std::size_t N,
    class=enable_if_t<N < typelist::size()>,
    class T
  > explicit variant (size<N>&&, ::std::true_type&&, T&& value) :
    data { }, tag { N }
  { ::new (this->target()) element<N> (::core::forward<T>(value)); }

  template <class T>
  using select_index = meta::either<
    meta::count<typelist, decay_t<T>>(),
    meta::index_of_t<typelist, decay_t<T>>,
    size<0>
  >;


public:

  /* The conditional_t used here allows us to first check if a given type
   * is declared in the variant and if it is, we will try to find its
   * constructor and immediately jump there, otherwise, we go the slower
   * route of trying to construct something from the value given.
   *
   * While this route is 'slower' this is a compile time performance issue and
   * will not impact runtime performance.
   *
   * Unfortunately we *do* instantiate templates several times, but there's
   * not much we can do about it.
   */
  template <
    class T,
    meta::inhibit<::std::is_same<decay_t<T>, variant>::value> = __LINE__
  > variant (T&& value) :
    variant {
      select_index<T> { },
      ::std::is_constructible<element<select_index<T>::value>, T> { },
      ::core::forward<T>(value)
    }
  { }

  template <
    class T,
    class... Args,
    meta::require<meta::count<typelist, T>() == 1> = __LINE__,
    meta::require<::std::is_constructible<T, Args...>::value> = __LINE__
  > variant (emplace_type_t<T>, Args&&... args) :
    data { },
    tag { meta::index_of<typelist, T>() }
  { ::new (this->target()) T(::core::forward<Args>(args)...); }

  template <
    ::std::size_t I,
    class... Args,
    meta::require<I < typelist::size()> = __LINE__,
    meta::require<
      ::std::is_constructible<element<I>, Args...>::value
    > = __LINE__
  > variant (emplace_index_t<I>, Args&&... args) :
    data { },
    tag { I }
  { ::new (this->target()) element<I>(::core::forward<Args>(args)...); }

  variant (variant const& that) :
    data { },
    tag { that.tag }
  { that.visit(copier { this->target() }); }

  variant (variant&& that) noexcept :
    data { }, tag { that.tag }
  { that.visit(mover { this->target() }); }

  variant () : variant { element<0> { } } { }

  ~variant () { this->visit(destroyer { }); }

  template <
    class T,
    meta::inhibit<::std::is_same<decay_t<T>, variant>::value> = __LINE__
  > variant& operator = (T&& value) {
    variant { ::core::forward<T>(value) }.swap(*this);
    return *this;
  }

  variant& operator = (variant const& that) {
    variant { that }.swap(*this);
    return *this;
  }

  variant& operator = (variant&& that) noexcept {
    this->visit(destroyer { });
    this->tag = that.tag;
    that.visit(mover { this->target() });
    return *this;
  }

  /* Placing these inside of the variant results in no implicit conversions
   * occuring with any potential constructor types.
   */
  bool operator == (variant const& that) const noexcept {
    if (this->tag != that.tag) { return false; }
    return that.visit(equality { this->target() });
  }

  bool operator < (variant const& that) const noexcept {
    if (this->tag != that.tag) { return this->tag < that.tag; }
    return that.visit(less_than { this->target() });
  }

  void swap (variant& that) noexcept(
    meta::all_of<typelist, is_nothrow_swappable>()
  ) {
    if (this->index() == that.index()) {
      that.visit(swapper { this->target() });
      return;
    }
    variant temp { ::core::move(*this) };
    *this = ::core::move(that);
    that = ::core::move(temp);
  }

  /* V stands for Visitor */
  template <class V, class... Args>
  auto visit (V&& visitor, Args&&... args) -> common_type_t<
    result_of_t<V(add_lvalue_reference_t<Ts>, Args...)>...
  > {
    using return_type = impl::result_t<
      V,
      meta::list<add_lvalue_reference_t<Ts>...>,
      Args...
    >;
    using function = add_pointer_t<return_type(V&&, void*, Args&&...)>;
    static auto const callers = make_array(
      impl::gen<V, Ts, void, function, Args...>()...
    );

    return callers[this->tag](
      ::core::forward<V>(visitor),
      this->target(),
      ::core::forward<Args>(args)...
    );
  }

  /* V stands for Visitor */
  template <class V, class... Args>
  auto visit (V&& visitor, Args&&... args) const -> common_type_t<
    result_of_t<V(add_lvalue_reference_t<add_const_t<Ts>>, Args...)>...
  > {
    using return_type = impl::result_t<
      V, meta::list<add_lvalue_reference_t<add_const_t<Ts>>...>, Args...
    >;

    using function = add_pointer_t<return_type(V&&, void const*, Args&&...)>;
    static auto const callers = make_array(
      impl::gen<V, add_const_t<Ts>, void const, function, Args...>()...
    );

    return callers[this->tag](
      ::core::forward<V>(visitor),
      this->target(),
      ::core::forward<Args>(args)...
    );
  }

  template <class... Vs>
  auto match (Vs&&... vs) -> decltype(
    this->visit(impl::make_overload(::core::forward<Vs>(vs)...))
  ) { return this->visit(impl::make_overload(::core::forward<Vs>(vs)...)); }

  template <class... Vs>
  auto match (Vs&&... vs) const -> decltype(
    this->visit(impl::make_overload(::core::forward<Vs>(vs)...))
  ) { return this->visit(impl::make_overload(::core::forward<Vs>(vs)...)); }

  /* These functions are undocumented and should not be used outside of core */
  template <::std::size_t N>
  add_pointer_t<add_const_t<element<N>>> cast () const noexcept {
    return static_cast<add_pointer_t<add_const_t<element<N>>>>(this->target());
  }

  template <::std::size_t N>
  add_pointer_t<element<N>> cast () noexcept {
    return static_cast<add_pointer_t<element<N>>>(this->target());
  }

  type_info const& type () const noexcept {
    return *this->visit(typeinfo { });
  }

  // Boost.Variant compatibility
  ::std::size_t which () const noexcept { return this->index(); }

  ::std::size_t index () const noexcept { return this->tag; }
  bool empty () const noexcept { return false; }

private:
  void const* target () const noexcept { return as_void(this->data); }
  void* target () noexcept { return as_void(this->data); }

  storage_type data;
  ::std::size_t tag;
};

template <class... Ts>
void swap (variant<Ts...>& lhs, variant<Ts...>& rhs) noexcept(
  noexcept(lhs.swap(rhs))
) { lhs.swap(rhs); }

template <::std::size_t I, class... Ts>
auto get (variant<Ts...> const* v) noexcept -> meta::when<
  I < sizeof...(Ts),
  add_pointer_t<add_const_t<meta::get<meta::list<Ts...>, I>>>
> {
  using t = add_pointer_t<add_const_t<meta::get<meta::list<Ts...>, I>>>;
  return v and v->index() == I
    ? static_cast<t>(v->template cast<I>())
    : nullptr;
}

template <::std::size_t I, class... Ts>
auto get (variant<Ts...>* v) noexcept -> meta::when<
  I < sizeof...(Ts),
  add_pointer_t<meta::get<meta::list<Ts...>, I>>
> {
  using t = add_pointer_t<meta::get<meta::list<Ts...>, I>>;
  return v and v->index() == I
    ? static_cast<t>(v->template cast<I>())
    : nullptr;
}

template <::std::size_t I, class... Ts>
auto get (variant<Ts...> const& v) noexcept(false) -> meta::when<
  I < sizeof...(Ts),
  add_lvalue_reference_t<add_const_t<meta::get<meta::list<Ts...>, I>>>
> {
  if (auto ptr = get<I>(::std::addressof(v))) { return *ptr; }
  throw_bad_variant_get();
}

template <::std::size_t I, class... Ts>
auto get (variant<Ts...>&& v) noexcept(false) -> meta::when<
  I < sizeof...(Ts),
  add_rvalue_reference_t<meta::get<meta::list<Ts...>, I>>
> {
  if (auto p = get<I>(::std::addressof(v))) { return ::core::move(*p); }
  throw_bad_variant_get();
}

template <::std::size_t I, class... Ts>
auto get (variant<Ts...>& v) noexcept(false) -> meta::when<
  I < sizeof...(Ts),
  add_lvalue_reference_t<meta::get<meta::list<Ts...>, I>>
> {
  if (auto ptr = get<I>(::std::addressof(v))) { return *ptr; }
  throw_bad_variant_get();
}

template <class T, class... Ts>
auto get (variant<Ts...> const* v) noexcept -> meta::unless<
  meta::find<meta::list<Ts...>, T>::empty(),
  decltype(get<meta::index_of<meta::list<Ts...>, T>()>(v))
> { return get<meta::index_of<meta::list<Ts...>, T>()>(v); }

template <class T, class... Ts>
auto get (variant<Ts...>* v) noexcept -> meta::unless<
  meta::find<meta::list<Ts...>, T>::empty(),
  decltype(get<meta::index_of<meta::list<Ts...>, T>()>(v))
> { return get<meta::index_of<meta::list<Ts...>, T>()>(v); }

template <class T, class... Ts>
auto get (variant<Ts...> const& v) noexcept(false) -> meta::unless<
  meta::find<meta::list<Ts...>, T>::empty(),
  decltype(get<meta::index_of<meta::list<Ts...>, T>()>(v))
> { return get<meta::index_of<meta::list<Ts...>, T>()>(v); }

template <class T, class... Ts>
auto get (variant<Ts...>&& v) noexcept(false) -> meta::unless<
  meta::find<meta::list<Ts...>, T>::empty(),
  decltype(core::move(get<meta::index_of<meta::list<Ts...>, T>()>(v)))
> { return core::move(get<meta::index_of<meta::list<Ts...>, T>()>(v)); }

template <class T, class... Ts>
auto get (variant<Ts...>& v) noexcept(false) -> meta::unless<
  meta::find<meta::list<Ts...>, T>::empty(),
  decltype(get<meta::index_of<meta::list<Ts...>, T>()>(v))
> { return get<meta::index_of<meta::list<Ts...>, T>()>(v); }

}} /* namespace core::v2 */

namespace std {

template <class... Ts>
struct hash<::core::v2::variant<Ts...>> {
  size_t operator () (::core::v2::variant<Ts...> const& value) const {
    return value.match(hash<Ts> { }...);
  }
};

} /* namespace std */

#endif /* CORE_VARIANT_HPP */
