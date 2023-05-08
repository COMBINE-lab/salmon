#ifndef SAFE_TYPE_SAFE_TYPE_HPP
#define SAFE_TYPE_SAFE_TYPE_HPP

#include <functional>
#include <istream>
#include <ostream>
#include <type_traits>
#include <utility>

#if __cplusplus >= 201703L
#define STRONG_NODISCARD [[nodiscard]]
#else
#define STRONG_NODISCARD
#endif

namespace strong
{

namespace impl
{
  template <typename T, typename ... V>
  using WhenConstructible = std::enable_if_t<std::is_constructible<T, V...>::value>;
}

template <typename M, typename T>
using modifier = typename M::template modifier<T>;

struct default_constructible
{
  template <typename T>
  class modifier
  {
  };
};

namespace impl {
  template <typename T>
  constexpr bool supports_default_construction(const ::strong::default_constructible::modifier<T>*)
  {
    return true;
  }
}

template <typename T, typename Tag, typename ... M>
class type : public modifier<M, type<T, Tag, M...>>...
{
public:
  template <typename type_ = type,
            bool = impl::supports_default_construction(static_cast<type_*>(nullptr))>
  constexpr
  type()
    noexcept(noexcept(T{}))
  : val{}
  {
  }

  template <typename U,
    typename = impl::WhenConstructible<T, std::initializer_list<U>>>
  constexpr
  explicit
  type(
    std::initializer_list<U> us
  )
    noexcept(noexcept(T{us}))
  : val{us}
  {
  }
  template <typename ... U,
            typename = std::enable_if_t<std::is_constructible<T, U&&...>::value && (sizeof...(U) > 0)>>
  constexpr
  explicit
  type(
    U&& ... u)
  noexcept(std::is_nothrow_constructible<T, U...>::value)
  : val(std::forward<U>(u)...)
  {}

  friend void swap(type& a, type& b) noexcept(
                                        std::is_nothrow_move_constructible<type>::value &&
                                        std::is_nothrow_move_assignable<type>::value
                                      )
  {
    using std::swap;
    swap(a.val, b.val);
  }

  STRONG_NODISCARD
  constexpr T& value_of() & noexcept { return val;}
  STRONG_NODISCARD
  constexpr const T& value_of() const & noexcept { return val;}
  STRONG_NODISCARD
  constexpr T&& value_of() && noexcept { return std::move(val);}
private:
  T val;
};

template <typename T>
struct is_safe_type : std::false_type {};
template <typename T, typename Tag, typename ... M>
struct is_safe_type<type<T, Tag, M...>> : std::true_type {};

namespace impl {
  template <typename T>
  using WhenSafeType = std::enable_if_t<is_safe_type<std::decay_t<T>>::value>;
  template <typename T>
  using WhenNotSafeType = std::enable_if_t<!is_safe_type<std::decay_t<T>>::value>;
}

template <typename T>
struct underlying_type
{
  using type = T;
};

template <typename T, typename Tag, typename ... M>
struct underlying_type<type<T, Tag, M...>>
{
  using type = T;
};

template <typename T>
using underlying_type_t = typename underlying_type<T>::type;

template <
  typename T,
  typename = impl::WhenSafeType<T>>
STRONG_NODISCARD
constexpr
auto
value_of(T&& t)
noexcept
-> decltype(std::forward<T>(t).value_of())
{
  return std::forward<T>(t).value_of();
}

template <
  typename T,
  typename = impl::WhenNotSafeType<T>>
constexpr
T&&
value_of(T&& t)
noexcept
{
  return std::forward<T>(t);
}

namespace impl
{
template <typename T, typename U>
static constexpr T &get(U &u) noexcept
{
  static_assert(is_safe_type<T>::value, "");
  static_assert(std::is_base_of<U, T>::value, "");
  return static_cast<T &>(u);
}

template <typename T, typename U>
inline constexpr const T &get(const U &u) noexcept
{
  static_assert(is_safe_type<T>::value, "");
  static_assert(std::is_base_of<U, T>::value, "");
  return static_cast<const T &>(u);
}

template <typename T, typename U>
inline constexpr T &&get(U &&u) noexcept
{
  static_assert(is_safe_type<T>::value, "");
  static_assert(std::is_base_of<std::remove_reference_t<U>, T>::value, "");
  return static_cast<T &&>(static_cast<U&&>(u));
}

template <typename T, typename U>
inline constexpr decltype(auto) access(U &&t) noexcept
{
  static_assert(is_safe_type<T>::value, "");
  static_assert(std::is_base_of<std::remove_reference_t<U>, T>::value, "");
  return get<T>(std::forward<U>(t)).value_of();
}

template <typename T>
struct type_of;

template <typename T, typename Tag, typename ... M>
struct type_of<::strong::type<T, Tag, M...>>
{
  using type = T;
};

template <typename T, typename Tag, typename ... M>
struct type_of<::strong::type<T, Tag, M...> &>
{
  using type = T &;
};

template <typename T, typename Tag, typename ... M>
struct type_of<::strong::type<T, Tag, M...> &&>
{
  using type = T &&;
};

template <typename T, typename Tag, typename ... M>
struct type_of<const ::strong::type<T, Tag, M...> &>
{
  using type = const T &;
};

template <typename T>
typename type_of<T>::type type_of_t();
}

struct equality
{
  template <typename T>
  class modifier;
};


template <typename T, typename Tag, typename ... M>
class equality::modifier<::strong::type<T, Tag, M...>>
{
  using type = ::strong::type<T, Tag, M...>;
public:
  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator==(
    const modifier& lh,
    const modifier& rh)
  noexcept(noexcept(std::declval<const T&>() == std::declval<const T&>()))
  -> decltype(std::declval<const T&>() == std::declval<const T&>())
  {
    return impl::access<type>(lh) == impl::access<type>(rh);
  }

  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator!=(
    const modifier& lh,
    const modifier& rh)
  noexcept(noexcept(std::declval<const T&>() != std::declval<const T&>()))
  -> decltype(std::declval<const T&>() != std::declval<const T&>())
  {
    return impl::access<type>(lh) != impl::access<type>(rh);
  }
};

namespace impl
{
  template <typename T>
  struct require_copy_constructible
  {
    static constexpr bool value = std::is_copy_constructible<underlying_type_t<T>>::value;
    static_assert(value, "underlying type must be copy constructible");
  };
  template <typename T>
  struct require_move_constructible
  {
    static constexpr bool value = std::is_move_constructible<underlying_type_t<T>>::value;
    static_assert(value, "underlying type must be move constructible");
  };
  template <typename T>
  struct require_copy_assignable
  {
    static constexpr bool value = std::is_copy_assignable<underlying_type_t<T>>::value;
    static_assert(value, "underlying type must be copy assignable");
  };
  template <typename T>
  struct require_move_assignable
  {
    static constexpr bool value = std::is_move_assignable<underlying_type_t<T>>::value;
    static_assert(value, "underlying type must be move assignable");
  };

  template <bool> struct valid_type;
  template <>
  struct valid_type<true> {};

  template <typename T>
  struct require_semiregular
    : valid_type<require_copy_constructible<T>::value &&
                 require_move_constructible<T>::value &&
                 require_copy_assignable<T>::value &&
                 require_move_assignable<T>::value>
  {
  };

}
struct semiregular
{
  template <typename>
  class modifier;
};

template <typename T, typename Tag, typename ... M>
class semiregular::modifier<::strong::type<T, Tag, M...>>
  : public default_constructible::modifier<T>
  , private impl::require_semiregular<T>
{
};

struct regular
{
  template <typename T>
  class modifier
    : public semiregular::modifier<T>
    , public equality::modifier<T>
  {
  };
};

struct ordered
{
  template <typename T>
  class modifier;
};


template <typename T, typename Tag, typename ... M>
class ordered::modifier<::strong::type<T, Tag, M...>>
{
  using type = ::strong::type<T, Tag, M...>;
public:
  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator<(
    const modifier& lh,
    const modifier& rh)
  noexcept(noexcept(std::declval<const T&>() < std::declval<const T&>()))
  -> decltype(std::declval<const T&>() < std::declval<const T&>())
  {
    return impl::access<type>(lh) < impl::access<type>(rh);
  }

  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator<=(
    const modifier& lh,
    const modifier& rh)
  noexcept(noexcept(std::declval<const T&>() <= std::declval<const T&>()))
  -> decltype(std::declval<const T&>() <= std::declval<const T&>())
  {
    return impl::access<type>(lh) <= impl::access<type>(rh);
  }

  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator>(
    const modifier& lh,
    const modifier& rh)
  noexcept(noexcept(std::declval<const T&>() > std::declval<const T&>()))
  -> decltype(std::declval<const T&>() > std::declval<const T&>())
  {
    return impl::access<type>(lh) > impl::access<type>(rh);
  }

  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator>=(
    const modifier& lh,
    const modifier& rh)
  noexcept(noexcept(std::declval<const T&>() >= std::declval<const T&>()))
  -> decltype(std::declval<const T&>() >= std::declval<const T&>())
  {
    return impl::access<type>(lh) >= impl::access<type>(rh);
  }
};

struct ostreamable
{
  template <typename T>
  class modifier
  {
  public:
    friend
    std::ostream&
    operator<<(
      std::ostream &os,
      const modifier &t)
    {
      return os << impl::access<T>(t);
    }
  };
};

struct istreamable
{
  template <typename T>
  class modifier
  {
  public:
    friend
    std::istream&
    operator>>(
      std::istream &is,
      modifier &t)
    {
      return is >> impl::access<T>(t);
    }
  };
};

struct iostreamable
{
  template <typename T>
  class modifier
    : public ostreamable::modifier<T>
    , public istreamable::modifier<T>
  {
  };
};

struct incrementable
{
  template <typename T>
  class modifier
  {
  public:
    constexpr
    T&
    operator++()
    noexcept(noexcept(++std::declval<T&>().value_of()))
    {
      auto &self = impl::get<T>(*this);
      ++value_of(self);
      return self;
    }

    constexpr
    T
    operator++(int)
    {
      T rv{impl::get<T>(*this)};
      ++*this;
      return rv;
    }
  };
};

struct decrementable
{
  template <typename T>
  class modifier
  {
  public:
    constexpr
    T&
    operator--()
    noexcept(noexcept(--std::declval<T&>().value_of()))
    {
      auto &self = impl::get<T>(*this);
      --value_of(self);
      return self;
    }

    constexpr
    T
    operator--(int)
    {
      T rv{impl::get<T>(*this)};
      --*this;
      return rv;
    }
  };
};

struct bicrementable
{
  template <typename T>
  class modifier
    : public incrementable::modifier<T>
    , public decrementable::modifier<T>
  {
  };
};

struct boolean
{
  template <typename T>
  class modifier
  {
  public:
    explicit constexpr operator bool() const
    noexcept(noexcept(static_cast<bool>(impl::type_of_t<T>())))
    {
      return static_cast<bool>(impl::access<T>(*this));
    }
  };
};

struct hashable
{
  template <typename T>
  class modifier{};
};

struct difference
{
  template <typename T>
  class modifier;
};

template <typename T, typename Tag, typename ... M>
class difference::modifier<::strong::type<T, Tag, M...>>
: public ordered::modifier<::strong::type<T, Tag, M...>>
{
  using type = ::strong::type<T, Tag, M...>;
public:
  type& operator+=(const type& t)
  noexcept(noexcept(std::declval<T&>() += value_of(t)))
  {
    impl::access<type>(*this) += value_of(t);
    return impl::get<type>(*this);
  }

  type& operator-=(const type& t)
    noexcept(noexcept(std::declval<T&>() -= value_of(t)))
  {
    impl::access<type>(*this) -= value_of(t);
    return impl::get<type>(*this);
  }

  type& operator*=(const T& t)
  noexcept(noexcept(std::declval<T&>() *= t))
  {
    impl::access<type>(*this) *= t;
    return impl::get<type>(*this);
  }

  type& operator/=(const T& t)
    noexcept(noexcept(std::declval<T&>() /= t))
  {
    impl::access<type>(*this) /= t;
    return impl::get<type>(*this);
  }

  friend
  type operator+(type lh, const type& rh)
  {
    lh += rh;
    return lh;
  }

  friend
  type operator-(type lh, const type& rh)
  {
    lh -= rh;
    return lh;
  }

  friend
  type operator*(type lh, const T& rh)
  {
    lh *= rh;
    return lh;
  }

  friend
  type operator*(const T& lh, type rh)
  {
    rh *= lh;
    return rh;
  }

  friend
  type operator/(type lh, const T& rh)
  {
    lh /= rh;
    return lh;
  }

  friend
  T operator/(const type& lh, const type& rh)
  {
    return value_of(lh) / value_of(rh);
  }
};

template <typename D>
struct affine_point
{
  template <typename T>
  class modifier;
};

namespace impl
{
  template <typename ...>
  using void_t = void;

  template <typename T, typename = void>
  struct subtractable : std::false_type {};

  template <typename T>
  struct subtractable<T, void_t<decltype(std::declval<const T&>() - std::declval<const T&>())>>
  : std::true_type {};
}


template <typename D>
template <typename T, typename Tag, typename ... M>
class affine_point<D>::modifier<::strong::type<T, Tag, M...>>
{
  using type = ::strong::type<T, Tag, M...>;
  static_assert(impl::subtractable<T>::value, "it must be possible to subtract instances of your underlying type");
  using diff_type = decltype(std::declval<const T&>() - std::declval<const T&>());
  static_assert(std::is_constructible<D, diff_type>::value,"");
public:
  STRONG_NODISCARD
  friend
  constexpr
  D
  operator-(
    const modifier& lh,
    const modifier& rh)
  {
    return D(impl::access<type>(lh) - impl::access<type>(rh));
  }

  type&
  operator+=(
    const D& d)
  noexcept(noexcept(std::declval<T&>() += value_of(d)))
  {
    impl::access<type>(*this) += value_of(d);
    return impl::get<type>(*this);
  }

  type&
  operator-=(
    const D& d)
  noexcept(noexcept(std::declval<T&>() -= value_of(d)))
  {
    impl::access<type>(*this) -= value_of(d);
    return impl::get<type>(*this);
  }

  STRONG_NODISCARD
  friend
  type
  operator+(
    const modifier& lh,
    const D& d)
  {
    return type(impl::access<type>(lh) + value_of(d));
  }

  STRONG_NODISCARD
  friend
  type
  operator+(
    const D& d,
    const modifier& rh)
  {
    return type(value_of(d) + impl::access<type>(rh));
  }

  STRONG_NODISCARD
  friend
  type
  operator-(
    const modifier& lh,
    const D& d)
  {
    return type(impl::access<type>(lh) - value_of(d));
  }
};

struct pointer
{
  template <typename T>
  class modifier;
};

template <typename T, typename Tag, typename ... M>
class pointer::modifier<::strong::type<T, Tag, M...>>
{
  using type = strong::type<T, Tag, M...>;
public:
  template <typename U = type, typename TT = T>
  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator==(
    const type& t,
    std::nullptr_t)
  noexcept(noexcept(std::declval<const TT&>() == nullptr))
  -> decltype(std::declval<const TT&>() == nullptr)
  {
    return value_of(t) == nullptr;
  }

  template <typename U = type, typename TT = T>
  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator==(
    std::nullptr_t,
    const type& t)
  noexcept(noexcept(nullptr == std::declval<const TT&>()))
  -> decltype(nullptr == std::declval<const TT&>())
  {
    return value_of(t) == nullptr;
  }

  template <typename U = type, typename TT = T>
  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator!=(
    const type& t,
    std::nullptr_t)
  noexcept(noexcept(std::declval<const TT&>() != nullptr))
  -> decltype(std::declval<const TT&>() != nullptr)
  {
    return value_of(t) != nullptr;
  }

  template <typename U = type, typename TT = T>
  STRONG_NODISCARD
  friend
  constexpr
  auto
  operator!=(
    std::nullptr_t,
    const type& t)
  noexcept(noexcept(nullptr != std::declval<const TT&>()))
  -> decltype(nullptr != std::declval<const TT&>())
  {
    return value_of(t) != nullptr;
  }

  STRONG_NODISCARD
  constexpr
  decltype(auto)
  operator*()
  const
  {
    return *impl::access<type>(*this);
  }

  constexpr decltype(auto) operator->() const { return &operator*();}
};

struct arithmetic
{
  template <typename T>
  class modifier
  {
  public:
    STRONG_NODISCARD
    friend
    constexpr
    T
    operator-(
      const T &lh)
    {
      return T{-value_of(lh)};
    }

    friend
    constexpr
    T&
    operator+=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) += value_of(rh)))
    {
      value_of(lh) += value_of(rh);
      return lh;
    }

    friend
    constexpr
    T&
    operator-=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) -= value_of(rh)))
    {
      value_of(lh) -= value_of(rh);
      return lh;
    }

    friend
    constexpr
    T&
    operator*=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) *= value_of(rh)))
    {
      value_of(lh) *= value_of(rh);
      return lh;
    }

    friend
    constexpr
    T&
    operator/=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) /= value_of(rh)))
    {
      value_of(lh) /= value_of(rh);
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator+(
      T lh,
      const T &rh)
    {
      lh += rh;
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator-(
      T lh,
      const T &rh)
    {
      lh -= rh;
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator*(
      T lh,
      const T &rh)
    {
      lh *= rh;
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator/(
      T lh,
      const T &rh)
    {
      lh /= rh;
      return lh;
    }
  };
};


struct bitarithmetic
{
  template <typename T>
  class modifier
  {
  public:
    friend
    constexpr
    T&
    operator&=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) &= value_of(rh)))
    {
      value_of(lh) &= value_of(rh);
      return lh;
    }

    friend
    constexpr
    T&
    operator|=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) |= value_of(rh)))
    {
      value_of(lh) |= value_of(rh);
      return lh;
    }

    friend
    constexpr
    T&
    operator^=(
      T &lh,
      const T &rh)
    noexcept(noexcept(value_of(lh) ^= value_of(rh)))
    {
      value_of(lh) ^= value_of(rh);
      return lh;
    }

    template <typename C>
    friend
    constexpr
    T&
    operator<<=(
      T &lh,
      C c)
    noexcept(noexcept(value_of(lh) <<= c))
    {
      value_of(lh) <<= c;
      return lh;
    }

    template <typename C>
    friend
    constexpr
    T&
    operator>>=(
      T &lh,
      C c)
    noexcept(noexcept(value_of(lh) >>= c))
    {
      value_of(lh) >>= c;
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator~(
      const T &lh)
    {
      auto v = value_of(lh);
      v = ~v;
      return T(v);
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator&(
      T lh,
      const T &rh)
    {
      lh &= rh;
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator|(
      T lh,
      const T &rh)
    {
      lh |= rh;
      return lh;
    }

    STRONG_NODISCARD
    friend
    constexpr
    T
    operator^(
      T lh,
      const T &rh)
    {
      lh ^= rh;
      return lh;
    }

    template <typename C>
    STRONG_NODISCARD
    friend
    constexpr
    T
    operator<<(
      T lh,
      C c)
    {
      lh <<= c;
      return lh;
    }

    template <typename C>
    STRONG_NODISCARD
    friend
    constexpr
    T
    operator>>(
      T lh,
      C c)
    {
      lh >>= c;
      return lh;
    }
  };
};
template <typename I = void>
struct indexed
{
  template <typename T>
  class modifier;
};


template <>
template <typename Type>
class indexed<void>::modifier
{
  using ref = typename impl::type_of<Type&>::type;
  using cref = typename impl::type_of<const Type&>::type;
  using rref = typename impl::type_of<Type&&>::type;
public:
  template <typename I>
  STRONG_NODISCARD
  auto
  operator[](
    const I& i)
  const &
  noexcept(noexcept(std::declval<cref>()[strong::value_of(i)]))
  -> decltype(std::declval<cref>()[strong::value_of(i)])
  {
    return impl::access<Type>(*this)[strong::value_of(i)];
  }

  template <typename I>
  STRONG_NODISCARD
  auto
  operator[](
    const I& i)
  &
  noexcept(noexcept(std::declval<ref>()[strong::value_of(i)]))
  -> decltype(std::declval<ref>()[strong::value_of(i)])
  {
    return impl::access<Type>(*this)[strong::value_of(i)];
  }

  template <typename I>
  STRONG_NODISCARD
  auto
  operator[](
    const I& i)
  &&
  noexcept(noexcept(std::declval<rref>()[strong::value_of(i)]))
  -> decltype(std::declval<rref>()[strong::value_of(i)])
  {
    return impl::access<Type>(*this)[strong::value_of(i)];
  }
  template <typename I, typename C = cref>
  STRONG_NODISCARD
  auto
  at(
    const I& i)
  const &
  -> decltype(std::declval<C>().at(strong::value_of(i)))
  {
    return impl::access<Type>(*this).at(strong::value_of(i));
  }
  template <typename I, typename R = ref>
  STRONG_NODISCARD
  auto
  at(
    const I& i)
  &
  -> decltype(std::declval<R>().at(strong::value_of(i)))
  {
    return impl::access<Type>(*this).at(strong::value_of(i));
  }
  template <typename I, typename R = rref>
  STRONG_NODISCARD
  auto
  at(
    const I& i)
  &&
  -> decltype(std::declval<R>().at(strong::value_of(i)))
  {
    return impl::access<Type>(*this).at(strong::value_of(i));
  }
};

template <typename I>
template <typename T, typename Tag, typename ... M>
class indexed<I>::modifier<type<T, Tag, M...>>
{
  using type = ::strong::type<T, Tag, M...>;
public:
  STRONG_NODISCARD
  auto
  operator[](
    const I& i)
  const &
  noexcept(noexcept(std::declval<const T&>()[strong::value_of(i)]))
  -> decltype(std::declval<const T&>()[strong::value_of(i)])
  {
    return impl::access<type>(*this)[strong::value_of(i)];
  }

  STRONG_NODISCARD
  auto
  operator[](
    const I& i)
  &
  noexcept(noexcept(std::declval<T&>()[strong::value_of(i)]))
  -> decltype(std::declval<T&>()[strong::value_of(i)])
  {
    return impl::access<type>(*this)[strong::value_of(i)];
  }

  STRONG_NODISCARD
  auto
  operator[](
    const I& i)
  &&
  noexcept(noexcept(std::declval<T&&>()[strong::value_of(i)]))
  -> decltype(std::declval<T&&>()[strong::value_of(i)])
  {
    return impl::access<type>(*this)[strong::value_of(i)];
  }

  STRONG_NODISCARD
  auto
  at(
    const I& i)
  const &
  -> decltype(std::declval<const T&>().at(strong::value_of(i)))
  {
    return impl::access<type>(*this).at(strong::value_of(i));
  }

  STRONG_NODISCARD
  auto
  at(
    const I& i)
  &
  -> decltype(std::declval<T&>().at(strong::value_of(i)))
  {
    return impl::access<type>(*this).at(strong::value_of(i));
  }

  STRONG_NODISCARD
  auto
  at(
    const I& i)
  &&
  -> decltype(std::declval<T&&>().at(strong::value_of(i)))
  {
    return impl::access<type>(*this).at(strong::value_of(i));
  }
};

class iterator
{
public:
  template <typename I, typename category = typename std::iterator_traits<I>::iterator_category>
  class modifier
    : public pointer::modifier<I>
      , public equality::modifier<I>
      , public incrementable::modifier<I>
  {
  };

  template <typename I>
  class modifier<I, std::bidirectional_iterator_tag>
    : public modifier<I, std::forward_iterator_tag>
      , public decrementable::modifier<I>
  {
  };
  template <typename I>
  class modifier<I, std::random_access_iterator_tag>
    : public modifier<I, std::bidirectional_iterator_tag>
      , public affine_point<typename std::iterator_traits<I>::difference_type>::template modifier<I>
      , public indexed<>::modifier<I>
      , public ordered::modifier<I>
  {
  };
};

class range
{
public:
  template <typename R>
  class modifier;
};

template <typename T, typename Tag, typename ... M>
class range::modifier<type<T, Tag, M...>>
{
  using type = ::strong::type<T, Tag, M...>;
  using r_iterator = decltype(std::declval<T&>().begin());
  using r_const_iterator = decltype(std::declval<const T&>().begin());
public:
  using iterator = ::strong::type<r_iterator, Tag, strong::iterator>;
  using const_iterator = ::strong::type<r_const_iterator, Tag, strong::iterator>;

  iterator
  begin()
  noexcept(noexcept(std::declval<T&>().begin()))
  {
    return iterator{impl::access<type>(*this).begin()};
  }

  iterator
  end()
  noexcept(noexcept(std::declval<T&>().end()))
  {
    return iterator{impl::access<type>(*this).end()};
  }

  const_iterator
  cbegin()
    const
  noexcept(noexcept(std::declval<const T&>().begin()))
  {
    return const_iterator{impl::access<type>(*this).begin()};
  }

  const_iterator
  cend()
    const
  noexcept(noexcept(std::declval<const T&>().end()))
  {
    return const_iterator{impl::access<type>(*this).end()};
  }

  const_iterator
  begin()
  const
  noexcept(noexcept(std::declval<const T&>().begin()))
  {
    return const_iterator{impl::access<type>(*this).begin()};
  }

  const_iterator
  end()
  const
  noexcept(noexcept(std::declval<const T&>().end()))
  {
    return const_iterator{impl::access<type>(*this).end()};
  }
};
}

namespace std {
template <typename T, typename Tag, typename ... M>
struct hash<::strong::type<T, Tag, M...>>
  : std::conditional_t<
    std::is_base_of<
      ::strong::hashable::modifier<
        ::strong::type<T, Tag, M...>
      >,
      ::strong::type<T, Tag, M...>
    >::value,
    hash<T>,
    std::false_type>
{
  using type = ::strong::type<T, Tag, M...>;
  decltype(auto)
  operator()(
    const ::strong::hashable::modifier<type>& t)
  const
  {
    return hash<T>::operator()(::strong::impl::access<type>(t));
  }
};
template <typename T, typename Tag, typename ... M>
struct is_arithmetic<::strong::type<T, Tag, M...>>
  : is_base_of<::strong::arithmetic::modifier<::strong::type<T, Tag, M...>>,
               ::strong::type<T, Tag, M...>>
{
};

template <typename T, typename Tag, typename ... M>
struct iterator_traits<::strong::type<T, Tag, M...>>
  : std::iterator_traits<T>
{
  using difference_type = typename std::iterator_traits<T>::difference_type;
  using value_type = typename std::iterator_traits<T>::value_type;
  using pointer = typename std::iterator_traits<T>::value_type;
  using reference = typename std::iterator_traits<T>::reference;
  using iterator_category = typename std::iterator_traits<T>::iterator_category;
};

}
#endif //SAFE_TYPE_SAFE_TYPE_HPP
