#ifndef CORE_MEMORY_RESOURCE_HPP
#define CORE_MEMORY_RESOURCE_HPP

#include <core/memory.hpp>

#include <atomic>
#include <new>

namespace core {
inline namespace v2 {
namespace pmr {

template <class T> struct polymorphic_allocator;

struct unsynchronized_pool_resource;
struct synchronized_pool_resource;
struct monotonic_buffer_resource;
struct memory_resource;

inline memory_resource* set_default_resource (memory_resource* mr) noexcept;
inline memory_resource* get_default_resource () noexcept;
inline memory_resource* null_memory_resource () noexcept;
inline memory_resource* new_delete_resource () noexcept;

inline bool operator == (
  memory_resource const&,
  memory_resource const&
) noexcept;

inline bool operator != (
  memory_resource const&,
  memory_resource const&
) noexcept;

template <class T>
bool operator == (
  polymorphic_allocator<T> const&,
  polymorphic_allocator<T> const&
) noexcept;

template <class T>
bool operator != (
  polymorphic_allocator<T> const&,
  polymorphic_allocator<T> const&
) noexcept;

}}} /* namespace core::v2::pmr */

namespace core {
inline namespace v2 {
namespace pmr {
namespace impl {

inline ::std::atomic<memory_resource*>& resource () noexcept {
  static ::std::atomic<memory_resource*> instance { new_delete_resource() };
  return instance;
}

template <class> struct resource_adaptor;

}}}} /* namespace core::v2::pmr::impl */

namespace core {
inline namespace v2 {
namespace pmr {

template <class T>
struct polymorphic_allocator {

  using value_type = T;

  template <class U>
  polymorphic_allocator (polymorphic_allocator<U> const& that) noexcept :
    mr { that.mr }
  { }

  polymorphic_allocator (memory_resource* mr) noexcept :
    mr { mr ? mr : get_default_resource() }
  { }

  polymorphic_allocator () noexcept :
    mr { get_default_resource() }
  { }

  memory_resource* resource () const noexcept { return this->mr; }

  value_type* allocate (::std::size_t n) {
    auto ptr = this->mr->allocate(n * sizeof(value_type), alignof(value_type));
    return static_cast<value_type*>(ptr);
  }

  void deallocate (value_type* ptr, ::std::size_t n) {
    this->mr->deallocate(ptr, n * sizeof(value_type), alignof(value_type));
  }

  polymorphic_allocator select_on_container_copy_construction () const {
    return polymorphic_allocator { };
  }

private:
  memory_resource* mr;
};

struct memory_resource {

  virtual ~memory_resource () noexcept = default;

  void* allocate (
    ::std::size_t bytes,
    ::std::size_t alignment = alignof(::std::max_align_t)
  ) noexcept(false) { return this->do_allocate(bytes, alignment); }

  void deallocate (
    void* ptr,
    ::std::size_t bytes,
    ::std::size_t alignment = alignof(::std::max_align_t)
  ) noexcept(false) { this->do_deallocate(ptr, bytes, alignment); }

  bool is_equal (memory_resource const& that) const noexcept {
    return this->do_is_equal(that);
  }

protected:

  virtual void do_deallocate (void*, ::std::size_t, ::std::size_t) = 0;
  virtual void* do_allocate (::std::size_t, ::std::size_t) = 0;
  virtual bool do_is_equal (memory_resource const&) const noexcept = 0;
};

struct monotonic_buffer_resource : memory_resource {

  virtual ~monotonic_buffer_resource () noexcept { }

  monotonic_buffer_resource& operator = (
    monotonic_buffer_resource const&
  ) = delete;

protected:

  virtual void do_deallocate (void*, ::std::size_t, ::std::size_t) override { }
};

inline memory_resource* set_default_resource (memory_resource* mr) noexcept {
  if (not mr) { mr = new_delete_resource(); }
  return impl::resource().exchange(mr);
}

inline memory_resource* get_default_resource () noexcept {
  return impl::resource();
}

inline memory_resource* null_memory_resource () noexcept {
  static struct : memory_resource {

    virtual void* do_allocate (::std::size_t, ::std::size_t) final {
      throw_bad_alloc();
    }

    virtual void do_deallocate (void*, ::std::size_t, ::std::size_t) final { }

    virtual bool do_is_equal (memory_resource const& mr) const noexcept final {
      return ::std::addressof(mr) == this;
    }

  } instance;
  return ::std::addressof(instance);
}

inline memory_resource* new_delete_resource () noexcept {
  static struct : memory_resource {

    virtual void* do_allocate (::std::size_t sz, ::std::size_t) final {
      return ::operator new (sz);
    }

    virtual void do_deallocate (void* p, ::std::size_t, ::std::size_t) final {
      return ::operator delete(p);
    }

    virtual bool do_is_equal (memory_resource const& mr) const noexcept final {
      return ::std::addressof(mr) == this;
    }

  } instance;
  return ::std::addressof(instance);
}

template <class Allocator>
using resource_adaptor = impl::resource_adaptor<
  typename ::std::allocator_traits<Allocator>::template rebind_alloc<char>
>;

inline bool operator == (
  memory_resource const& lhs,
  memory_resource const& rhs
) noexcept {
  return ::std::addressof(lhs) == ::std::addressof(rhs) or lhs.is_equal(rhs);
}

inline bool operator != (
  memory_resource const& lhs,
  memory_resource const& rhs
) noexcept { return not (lhs == rhs); }

template <class T>
bool operator == (
  polymorphic_allocator<T> const& lhs,
  polymorphic_allocator<T> const& rhs
) noexcept { return *lhs.resource() == *rhs.resource(); }

template <class T>
bool operator != (
  polymorphic_allocator<T> const& lhs,
  polymorphic_allocator<T> const& rhs
) noexcept { return *lhs.resource() != *rhs.resource(); }

}}} /* namespace core::v2::pmr */

namespace core {
inline namespace v2 {
namespace pmr {
namespace impl {

template <class Allocator>
struct resource_adaptor : memory_resource {

  using allocator_type = Allocator;
  using alloc_traits = ::std::allocator_traits<allocator_type>;

  static constexpr bool same_cp = ::std::is_same<
    typename alloc_traits::value_type const*,
    typename alloc_traits::const_pointer
  >::value;

  static constexpr bool same_p = ::std::is_same<
    typename alloc_traits::value_type*,
    typename alloc_traits::pointer
  >::value;

  static constexpr bool same_cv = ::std::is_same<
    typename alloc_traits::const_void_pointer,
    void const*
  >::value;

  static constexpr bool same_v = ::std::is_same<
    typename alloc_traits::void_pointer,
    void*
  >::value;

  static_assert(same_cp, "const_pointer must be same as value_type const*");
  static_assert(same_cv, "const_void_pointer must be same as void const*");
  static_assert(same_p, "pointer must be same as value_type*");
  static_assert(same_v, "void_pointer must be same as void*");

  template <class Alloc>
  resource_adaptor (Alloc&& alloc) :
    alloc { ::core::forward<Alloc>(alloc) }
  { }

  resource_adaptor (resource_adaptor const&) = default;
  resource_adaptor () = default;

  resource_adaptor& operator = (resource_adaptor const&) = default;

  allocator_type get_allocator () const { return this->alloc; }

private:

  virtual void* do_allocate (::std::size_t size, ::std::size_t) final {
    return alloc_traits::allocate(this->alloc, size);
  }

  virtual void do_deallocate (void* p, ::std::size_t sz, ::std::size_t) final {
    alloc_traits::deallocate(this->alloc, p, sz);
  }

  virtual bool do_is_equal (memory_resource const& that) const noexcept final {
    auto ptr = dynamic_cast<resource_adaptor const*>(&that);
    if (not ptr) { return false; }
    return ptr->alloc == this->alloc;
  }

  allocator_type alloc;
};

}}}} /* namespace core::v2::pmr::impl */

#endif /* CORE_MEMORY_RESOURCE_HPP */
