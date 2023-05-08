// itlib-pod-vector v1.04
//
// A vector of PODs. Similar to std::vector, but doesn't call constructors or
// destructors and instead uses memcpy and memmove to manage the data
//
// SPDX-License-Identifier: MIT
// MIT License:
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
//  1.04 (2021-08-05) Bugfix! Fixed return value of erase
//  1.03 (2021-06-08) Prevent memcmp calls with nullptr
//  1.02 (2021-06-08) Noexcept move ctor and move assignment operator
//  1.01 (2020-10-28) Switched static assert from is_pod to is_trivial
//  1.00 (2020-10-18) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::pod_vector, which similar to std::vector:
// * It keeps the data in a contiguous memory block
// * Has the same public methods and operators and features like random-access
// But:
// * Operates only ot PODs
// * Doesn't call constructors, destructors, move and assign operators
// * Instead uses memcpy and memmove to manage the data
// Thus, it achieves a much better performance, especially in Debug mode.
//
// pod_vector also allows "recast" where you can convert pod_vector<T> to
// pod_vector<U>. This is very useful when operating with signed/unsigned char
// for example.
//
// except for the methods which are the same as std::vector, itlib::pod_vector
// also provides the following:
// * size_t byte_size() const; - size of data in bytes
// * recast_copy_from(other_vec) - copies from other vec. Note that this will
//   lose data if the byte size of other_vec's data is not divisible by
//   sizeof(T)
// * recast_take_from(other_vec) - moves from other vec. Note that this will
//   lose data if the byte size of other_vec's data is not divisible by
//   sizeof(T)
//
// pod_vector uses pod_allocator, which needs to have methods to allocate,
// deallocate, and reallocate. The default version uses malloc, free, and
// realloc. If you make your own allocator you must conform to the definitons
// of these functions.
// The allocator must provide the following interface:
// * using size_type = ...; - size type for allocator and vector
// * void* malloc(size_type size); - allocate memory
// * void* realloc(void* old, size_type new_size); - allocate/reallocate memory
// * void free(void* mem); - free memory which was allocated here
// * size_type max_size(); - max available memory
// * bool zero_fill_new(); - whether pod_vector should to zerofill new elements
// * size_type realloc_wasteful_copy_size() - when to use reallocate on grows
//
//                  TESTS
//
// You can find unit tests for pod_vector in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <type_traits>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <cstring>

namespace itlib
{

namespace impl
{
class pod_allocator
{
public:
    using size_type = size_t;
    static void* malloc(size_type size) { return std::malloc(size); }
    static void* realloc(void* old, size_type new_size) { return std::realloc(old, new_size); }
    static void free(void* mem) { std::free(mem); }
    static constexpr size_type max_size() { return ~size_type(0); }
    static constexpr bool zero_fill_new() { return true; }
    static constexpr size_type realloc_wasteful_copy_size() { return 4096; }
};
}

template<typename T, class Alloc = impl::pod_allocator>
class pod_vector
{
    static_assert(std::is_trivial<T>::value, "itlib::pod_vector with non-trivial type");

    template<typename U, typename A>
    friend class pod_vector; // so we can move between types

public:
    using allocator_type = Alloc;
    using value_type = T;
    using size_type = typename Alloc::size_type;
    using reference = T&;
    using const_reference = const T&;
    using pointer = T*;
    using const_pointer = const T*;
    using iterator = pointer;
    using const_iterator = const_pointer;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    pod_vector()
        : pod_vector(Alloc())
    {}

    explicit pod_vector(const Alloc& alloc)
        : m_capacity(0)
        , m_alloc(alloc)
    {
        m_begin = m_end = nullptr;
    }

    explicit pod_vector(size_t count, const Alloc& alloc = Alloc())
        : pod_vector(alloc)
    {
        resize(count);
    }

    pod_vector(size_t count, const T& value, const Alloc& alloc = Alloc())
        : pod_vector(alloc)
    {
        assign_fill(count, value);
    }

    template <typename InputIterator, typename = decltype(*std::declval<InputIterator>())>
    pod_vector(InputIterator first, InputIterator last, const Alloc& alloc = Alloc())
        : pod_vector(alloc)
    {
        assign_copy(first, last);
    }

    pod_vector(std::initializer_list<T> l, const Alloc& alloc = Alloc())
        : pod_vector(alloc)
    {
        assign_copy(l.begin(), l.end());
    }

    pod_vector(const pod_vector& other)
        : pod_vector(other, other.get_allocator())
    {}

    pod_vector(const pod_vector& other, const Alloc& alloc)
        : pod_vector(alloc)
    {
        assign_copy(other.begin(), other.end());
    }

    pod_vector(pod_vector&& other) noexcept
        : m_begin(other.m_begin)
        , m_end(other.m_end)
        , m_capacity(other.m_capacity)
        , m_alloc(std::move(other.m_alloc))
    {
        other.m_begin = other.m_end = nullptr;
        other.m_capacity = 0;
    }

    ~pod_vector()
    {
        m_alloc.free(m_begin);
    }

    pod_vector& operator=(const pod_vector& other)
    {
        if (this == &other) return *this; // prevent self usurp
        clear();
        assign_copy(other.begin(), other.end());
        return *this;
    }

    pod_vector& operator=(pod_vector&& other) noexcept
    {
        if (this == &other) return *this; // prevent self usurp

        m_alloc.free(m_begin);

        m_alloc = std::move(other.m_alloc);
        m_capacity = other.m_capacity;
        m_begin = other.m_begin;
        m_end = other.m_end;

        other.m_begin = other.m_end = nullptr;
        other.m_capacity = 0;

        return *this;
    }

    template <typename U, typename UAlloc>
    void recast_copy_from(const pod_vector<U, UAlloc>& other)
    {
        clear();
        auto new_size = other.byte_size() / sizeof(T);
        auto cast = reinterpret_cast<const T*>(other.data());
        assign_copy(cast, cast + new_size);
    }

    template <typename U, typename UAlloc>
    void recast_take_from(pod_vector<U, UAlloc>&& other)
    {
        m_alloc.free(m_begin);

        auto new_size = other.byte_size() / sizeof(T);
        auto cast = reinterpret_cast<T*>(other.data());
        m_begin = cast;
        m_end = m_begin + new_size;

        m_capacity = (sizeof(U) * other.capacity()) / sizeof(T);

        m_alloc = std::move(other.m_alloc);

        other.m_begin = other.m_end = nullptr;
        other.m_capacity = 0;
    }

    void assign(size_type count, const T& value)
    {
        assign_fill(count, value);
    }

    template <typename InputIterator, typename = decltype(*std::declval<InputIterator>())>
    void assign(InputIterator first, InputIterator last)
    {
        assign_copy(first, last);
    }

    void assign(std::initializer_list<T> ilist)
    {
        assign_copy(ilist.begin(), ilist.end());
    }

    allocator_type get_allocator() const noexcept
    {
        return m_alloc;
    }

    const_reference at(size_type i) const
    {
        return *(m_begin + i);
    }

    reference at(size_type i)
    {
        return *(m_begin + i);
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
        return *(m_end - 1);
    }

    reference back()
    {
        return *(m_end - 1);
    }

    const_pointer data() const noexcept
    {
        return m_begin;
    }

    pointer data() noexcept
    {
        return m_begin;
    }

    // iterators
    iterator begin() noexcept
    {
        return m_begin;
    }

    const_iterator begin() const noexcept
    {
        return m_begin;
    }

    const_iterator cbegin() const noexcept
    {
        return m_begin;
    }

    iterator end() noexcept
    {
        return m_end;
    }

    const_iterator end() const noexcept
    {
        return m_end;
    }

    const_iterator cend() const noexcept
    {
        return m_end;
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
        return m_begin == m_end;
    }

    size_type size() const noexcept
    {
        return m_end - m_begin;
    }

    size_type max_size() const noexcept
    {
        return m_alloc.max_size();
    }

    size_t byte_size() const noexcept
    {
        return sizeof(value_type) * size();
    }

    void reserve(size_t desired_capacity)
    {
        if (desired_capacity <= m_capacity) return;

        auto new_cap = get_new_capacity(desired_capacity);
        auto s = size();

        T* new_buf;
        if (m_capacity - s > m_alloc.realloc_wasteful_copy_size())
        {
            // we decided that it would be more wasteful to use realloc and
            // copy more than needed than it would be to malloc and manually copy
            new_buf = pointer(m_alloc.malloc(new_cap * sizeof(T)));
            if (s) memcpy(new_buf, m_begin, byte_size());
            m_alloc.free(m_begin);
        }
        else
        {
            new_buf = pointer(m_alloc.realloc(m_begin, sizeof(value_type) * new_cap));
        }

        m_begin = new_buf;
        m_end = new_buf + s;
        m_capacity = new_cap;
    }

    size_t capacity() const noexcept
    {
        return m_capacity;
    }

    void shrink_to_fit()
    {
        const auto s = size();

        if (s == m_capacity) return;

        if (s == 0)
        {
            m_alloc.free(m_begin);
            m_capacity = 0;
            m_begin = m_end = nullptr;
            return;
        }

        auto new_buf = pointer(m_alloc.realloc(m_begin, sizeof(value_type) * s));

        m_begin = new_buf;
        m_end = new_buf + s;
        m_capacity = s;
    }

    // modifiers
    void clear() noexcept
    {
        m_end = m_begin;
    }

    iterator insert(const_iterator position, const value_type& val)
    {
        auto pos = grow_at(position, 1);
        *pos = val;
        return pos;
    }

    iterator insert(const_iterator position, size_type count, const value_type& val)
    {
        auto pos = grow_at(position, count);
        fill(pos, count, val);
        return pos;
    }

    template <typename InputIterator, typename = decltype(*std::declval<InputIterator>())>
    iterator insert(const_iterator position, InputIterator first, InputIterator last)
    {
        auto pos = grow_at(position, last - first);
        copy_not_aliased(pos, first, last);
        return pos;
    }

    iterator insert(const_iterator position, std::initializer_list<T> ilist)
    {
        auto pos = grow_at(position, ilist.size());
        copy_not_aliased(pos, ilist.begin(), ilist.end());
        return pos;
    }

    // for compatibility
    iterator emplace(const_iterator position, const_reference& val)
    {
        return insert(position, val);
    }

    iterator erase(const_iterator position)
    {
        return shrink_at(position, 1);
    }

    iterator erase(const_iterator first, const_iterator last)
    {
        return shrink_at(first, last - first);
    }

    // for compatibility
    reference emplace_back()
    {
        reserve(size() + 1);
        ++m_end;
        return back();
    }

    reference push_back(const_reference val)
    {
        return emplace_back() = val;
    }

    // for compatibility
    reference emplace_back(const_reference val)
    {
        return push_back(val);
    }

    void pop_back()
    {
        shrink_at(m_end - 1, 1);
    }

    void resize(size_type n, const value_type& val)
    {
        reserve(n);
        fill(m_end, n, val);
        m_end = m_begin + n;
    }

    void resize(size_type n)
    {
        reserve(n);
        if (n > size() && m_alloc.zero_fill_new()) zero_fill(m_end, n - size());
        m_end = m_begin + n;
    }

    void swap(pod_vector& other)
    {
        auto tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

private:
    static void zero_fill(T* p, size_type s)
    {
        std::memset(p, 0, s * sizeof(T));
    }

    // fill count elements from p with value
    static void fill(T* p, size_type count, const T& value)
    {
        // we could look for optimizations here
        // we could use memset if sizeof(T) == 1
        // however all observed compilers ended up recognizing this
        // and using memset with optimizations in such a case
        // that's enough for us
        for (size_type i=0; i<count; ++i)
        {
            *p++ = value;
        }
    }

    template <typename InputIterator>
    static void copy_not_aliased(T* p, InputIterator begin, InputIterator end)
    {
        // much like above when InputIterator is const T* or bitwise equivalent to T*
        // (say int/unsinged)
        // compilers with optimizations end up recognizing the case and using memcpy
        for (auto i = begin; i!=end; ++i)
        {
            *p++ = *i;
        }
    }

    // still for extra help, we can provide this (alsto it will be faster in debug)
    static void copy_not_aliased(T* p, const T* begin, const T* end)
    {
        auto s = size_t(end - begin) * sizeof(T);
        if (s == 0) return;
        std::memcpy(p, begin, s);
    }

    // calculate a new capacity so that it's at least desired_capacity
    size_type get_new_capacity(size_type desired_capacity) const
    {
        if (m_capacity == 0)
        {
            return desired_capacity;
        }
        else
        {
            auto new_cap = m_capacity;

            while (new_cap < desired_capacity)
            {
                // grow by roughly 1.5
                new_cap *= 3;
                ++new_cap;
                new_cap /= 2;
            }

            return new_cap;
        }
    }

    // increase the size by splicing the elements in such a way that
    // a hole of uninitialized elements is left at position, with size num
    // returns the (potentially new) address of the hole
    T* grow_at(const T* cp, size_type num)
    {
        const auto s = size();
        auto offset = cp - m_begin;

        if (cp == m_end)
        {
            resize(s + num);
        }
        else
        {
            if (s + num > m_capacity)
            {
                auto new_cap = get_new_capacity(s + num);
                T* new_buf;
                if (m_capacity - offset > m_alloc.realloc_wasteful_copy_size())
                {
                    // we decided that it would be more wasteful to use realloc and
                    // copy more than needed than it would be to malloc and manually copy
                    new_buf = pointer(m_alloc.malloc(new_cap * sizeof(T)));
                    if (offset) memcpy(new_buf, m_begin, offset * sizeof(T));
                    if (m_alloc.zero_fill_new()) zero_fill(new_buf + offset, num);
                    memcpy(new_buf + offset + num, m_begin + offset, (s - offset) * sizeof(T));
                    m_alloc.free(m_begin);
                }
                else
                {
                    new_buf = pointer(m_alloc.realloc(m_begin, sizeof(value_type) * new_cap));
                    std::memmove(new_buf + offset + num, new_buf + offset, (s - offset) * sizeof(T));
                    if (m_alloc.zero_fill_new()) zero_fill(new_buf + offset, num);
                }
                m_begin = new_buf;
                m_capacity = new_cap;
            }
            else
            {
                std::memmove(m_begin + offset + num, m_begin + offset, (s - offset) * sizeof(T));
                if (m_alloc.zero_fill_new()) zero_fill(m_begin + offset, num);
            }
        }

        m_end = m_begin + s + num;
        return m_begin + offset;
    }

    // remove elements from cp to cp+num, shifting the rest
    // returns one after the removed
    T* shrink_at(const T* cp, size_type num)
    {
        const auto s = size();
        if (s == num)
        {
            clear();
            return m_end;
        }

        auto position = const_cast<T*>(cp);

        std::memmove(position, position + num, size_t(m_end - position - num) * sizeof(T));

        m_end -= num;

        return position;
    }

    // grows buffer only on empty vectors
    void safe_grow(size_t capacity)
    {
        if (capacity <= m_capacity) return;

        m_alloc.free(m_begin);

        m_capacity = get_new_capacity(capacity);
        m_begin = m_end = pointer(m_alloc.malloc(sizeof(value_type) * m_capacity));
    }

    // fill empty vector with given value
    void assign_fill(size_type count, const T& value)
    {
        safe_grow(count);
        fill(m_begin, count, value);
        m_end = m_begin + count;
    }

    // fill empty vector with values from [first;last)
    template <class InputIterator>
    void assign_copy(InputIterator first, InputIterator last)
    {
        auto count = last - first;
        safe_grow(count);
        copy_not_aliased(m_begin, first, last);
        m_end = m_begin + count;
    }

    pointer m_begin;
    pointer m_end;

    size_t m_capacity;

    Alloc m_alloc;
};

template<typename T, class Alloc>
bool operator==(const pod_vector<T, Alloc>& a, const pod_vector<T, Alloc>& b)
{
    if (a.size() != b.size()) return false;
    if (a.empty()) return true;
    return std::memcmp(a.data(), b.data(), a.byte_size()) == 0;
}

template<typename T, class Alloc>
bool operator!=(const pod_vector<T, Alloc>& a, const pod_vector<T, Alloc>& b)
{
    if (a.size() != b.size()) return true;
    if (a.empty()) return false;
    return std::memcmp(a.data(), b.data(), a.byte_size()) != 0;
}

}
