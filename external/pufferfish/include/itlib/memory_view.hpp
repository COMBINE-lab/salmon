// itlib-memory-view v1.01
//
// A view of a chunk of memory which makes it look as a std::vector sans
// the size modifying functions
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2016-2017 Chobolabs Inc.
// Copyright(c) 2020-2021 Borislav Stanimirov
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files(the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and / or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions :
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT.IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
//                  VERSION HISTORY
//
//  1.01 (2021-10-07) Added slicing functionalities
//  1.00 (2020-10-14) Rebranded release from chobo-memory-view
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// A memory view is a class which attaches to an chunk of memory and provides
// a view to it (optionally changing its type) with a std::vector-like
// interface, excluding the methods which change a vector's size.
//
// The library includes two classes, for viewing const and non-const memory:
// memory_view, and const_memory_view. To automatically generate
// the appropriate pointer, use `make_memory_view(ptr, size)`.
// where size is the number of elements your view will have
//
// `make_memory_view` also has an overload for std::vector like types which
// can be useful if you want a view at a slise of the vector
// you can use it like this: `make_memory_view<OptionalNewType>(vector, offset = 0)`
// where the offset is from the beginning of the vector.
//
// It is the responsibility of the user to make sure the pointer given to a
// view remains valid throughout the view's lifetime. All methods of the library
// assume this.
//
// Copying a memory view will copy the unrerlying pointer (and not its data),
// resulting in another view to the same memory.
//
// To copy the data from one memory view to another, you can use std::copy
//
// Example:
//
// int* buf = new int[10];
// auto view = make_memory_view(buf, 10);
// for (auto& i : view) i = 0; // fill buf with zeroes
//
// Reference:
//
// memory_view has the most common std::vector methods and operators
// const_memory_view has the most common std::vector const methods and operators
// Both exclude the methods and operators which would result in changing the
// a vector's size.
//
// Besides those, both classes define a ptr-like interface:
//
// void reset(void* ptr = nullptr, size_t size = 0) - reset the internal pointer
// void reset(size_t size) - change the size, without changing the pointer
// operator bool() const - returns whether the internal pointer is valid
// T* get() noexcept - return the internal pointer (same as data())
// slice(size_t offset, size_t count = ~0) - returns a new memory view which is a slice
//   from this [offset, offset+count). If the slice would go beyond the end of the view
//   the end of the view is used
// remove_prefix(size_t n) - moves the start by n elements. UB if n > size
// remove_suffix(size_t n) - moves the end by n elements towards start. UB if n > size
//
//
//                  Configuration
//
// itlib::memory_view a single configurable setting:
//
//                  Config bounds checks:
//
// By default bounds checks are made in debug mode (via an assert) when accessing
// elements (with `at` or `[]`). Iterators are not checked (yet...)
//
// To disable them, you can define ITLIB_MEMORY_VIEW_NO_DEBUG_BOUNDS_CHECK
// before including the header.
//
//
//                  TESTS
//
// You can find unit tests for memory_view in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <cstddef>
#include <iterator>

#if defined(ITLIB_MEMORY_VIEW_NO_DEBUG_BOUNDS_CHECK)
#   define I_ITLIB_MEMORY_VIEW_BOUNDS_CHECK(i)
#else
#   include <cassert>
#   define I_ITLIB_MEMORY_VIEW_BOUNDS_CHECK(i) assert((i) < this->size())
#endif

namespace itlib
{

template <typename T>
class memory_view
{
public:

    using value_type = T;
    using size_type = size_t;
    using difference_type = ptrdiff_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    memory_view() = default;

    memory_view(T* ptr, size_t size)
        : m_ptr(ptr)
        , m_size(size)
    {}

    memory_view(const memory_view&) = default;
    memory_view(memory_view&&) = default;

    memory_view& operator=(const memory_view&) = default;
    memory_view& operator=(memory_view&&) = default;

    void reset(T* ptr = nullptr, size_t size = 0)
    {
        m_ptr = ptr;
        m_size = size;
    }

    void reset(size_t size)
    {
        m_size = size;
    }

    explicit operator bool() const
    {
        return !!m_ptr;
    }

    pointer get() noexcept { return m_ptr; }
    const_pointer get() const noexcept { return m_ptr; }

    const_reference at(size_type i) const
    {
        I_ITLIB_MEMORY_VIEW_BOUNDS_CHECK(i);
        return *(m_ptr + i);
    }

    reference at(size_type i)
    {
        I_ITLIB_MEMORY_VIEW_BOUNDS_CHECK(i);
        return *(m_ptr + i);
    }

    const_reference operator[](size_type i) const
    {
        return at(i);
    }

    reference operator[](size_type i)
    {
        return at(i);
    }

    const_reference front() const
    {
        return at(0);
    }

    reference front()
    {
        return at(0);
    }

    const_reference back() const
    {
        return at(m_size - 1);
    }

    reference back()
    {
        return at(m_size - 1);
    }

    const_pointer data() const noexcept
    {
        return m_ptr;
    }

    pointer data() noexcept
    {
        return m_ptr;
    }

    // iterators
    iterator begin() noexcept
    {
        return data();
    }

    const_iterator begin() const noexcept
    {
        return data();
    }

    const_iterator cbegin() const noexcept
    {
        return data();
    }

    iterator end() noexcept
    {
        return data() + m_size;
    }

    const_iterator end() const noexcept
    {
        return data() + m_size;
    }

    const_iterator cend() const noexcept
    {
        return data() + m_size;
    }

    reverse_iterator rbegin() noexcept
    {
        return reverse_iterator(end());
    }

    const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    const_reverse_iterator crbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    reverse_iterator rend() noexcept
    {
        return reverse_iterator(begin());
    }

    const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    const_reverse_iterator crend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    // capacity
    bool empty() const noexcept
    {
        return m_size == 0;
    }

    size_t size() const noexcept
    {
        return m_size;
    }

    // slicing
    memory_view slice(size_t off, size_t count = size_t(-1)) const
    {
        if (off > m_size) return memory_view(m_ptr + m_size, 0);
        auto newSize = m_size - off;
        if (count > newSize) count = newSize;
        return memory_view(m_ptr + off, count);
    }

    void remove_prefix(size_t n)
    {
        m_ptr += n;
        m_size -= n;
    }

    void remove_suffix(size_t n)
    {
        m_size -= n;
    }

private:
    T* m_ptr = nullptr;
    size_t m_size = 0;
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
class const_memory_view
{
public:

    using value_type = T;
    using size_type = size_t;
    using difference_type = ptrdiff_t;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using iterator = T*;
    using const_iterator = const T*;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    const_memory_view() = default;

    const_memory_view(const T* ptr, size_t size)
        : m_ptr(ptr)
        , m_size(size)
    {}

    const_memory_view(const memory_view<T>& view)
    {
        m_ptr = view.data();
        m_size = view.size();
    }

    const_memory_view(const const_memory_view&) = default;
    const_memory_view(const_memory_view&&) = default;

    const_memory_view& operator=(const const_memory_view&) = default;
    const_memory_view& operator=(const_memory_view&&) = default;

    void reset(const T* ptr = nullptr, size_t size = 0)
    {
        m_ptr = ptr;
        m_size = size;
    }

    void reset(size_t size)
    {
        m_size = size;
    }

    explicit operator bool() const
    {
        return !!m_ptr;
    }

    const_pointer get() const noexcept { return m_ptr; }

    const_reference at(size_type i) const
    {
        I_ITLIB_MEMORY_VIEW_BOUNDS_CHECK(i);
        return *(m_ptr + i);
    }

    const_reference operator[](size_type i) const
    {
        return at(i);
    }

    const_reference front() const
    {
        return at(0);
    }

    const_reference back() const
    {
        return at(m_size - 1);
    }

    const_pointer data() const noexcept
    {
        return m_ptr;
    }

    // iterators

    const_iterator begin() const noexcept
    {
        return data();
    }

    const_iterator cbegin() const noexcept
    {
        return data();
    }

    const_iterator end() const noexcept
    {
        return data() + m_size;
    }

    const_iterator cend() const noexcept
    {
        return data() + m_size;
    }

    const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    const_reverse_iterator crbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    const_reverse_iterator crend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    // capacity
    bool empty() const noexcept
    {
        return m_size == 0;
    }

    size_t size() const noexcept
    {
        return m_size;
    }

    // slicing
    const_memory_view slice(size_t off, size_t count = size_t(-1)) const
    {
        if (off > m_size) return const_memory_view(m_ptr + m_size, 0);
        auto newSize = m_size - off;
        if (count > newSize) count = newSize;
        return const_memory_view(m_ptr + off, count);
    }

    void remove_prefix(size_t n)
    {
        m_ptr += n;
        m_size -= n;
    }

    void remove_suffix(size_t n)
    {
        m_size -= n;
    }

private:
    const T* m_ptr = nullptr;
    size_t m_size = 0;
};

///////////////////////////////////////////////////////////////////////////////

template <typename T>
memory_view<T> make_memory_view(T* ptr, size_t size)
{
    return memory_view<T>(ptr, size);
}

template <typename T>
const_memory_view<T> make_memory_view(const T* ptr, size_t size)
{
    return const_memory_view<T>(ptr, size);
}

template <typename T, size_t N>
const_memory_view<T> make_memory_view(const T(&ar)[N])
{
    return const_memory_view<T>(ar, N);
}

template <typename Container>
auto make_memory_view(Container& c, size_t offset = 0) -> memory_view<typename Container::value_type>
{
    return memory_view<typename Container::value_type>(c.data() + offset, c.size() - offset);
}

template <typename Container>
auto make_memory_view(const Container& c, size_t offset = 0) -> const_memory_view<typename Container::value_type>
{
    return const_memory_view<typename Container::value_type>(c.data() + offset, c.size() - offset);
}

}
