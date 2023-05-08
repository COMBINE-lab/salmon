// itlib-static-vector v1.04
//
// std::vector-like class with a fixed capacity
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2016-2019 Chobolabs Inc.
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
//  1.04 (2021-11-18) Added assign ops
//  1.03 (2021-10-05) Don't rely on operator!= from T. Use operator== instead
//  1.02 (2021-08-04) emplace_back() returns a reference as per C++17
//  1.01 (2021-08-04) capacity() and max_size() to static constexpr methods
//  1.00 (2020-10-14) Rebranded release from chobo-static-vector
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::static_vector, which is almost a drop-in
// replacement of std::vector, but has a fixed capacity as a template argument.
// It gives you the benefits of using std::array (cache-locality) with the
// flexibility of std::vector - dynamic size.
// However the size may never exceed its initially set capacity.
//
// Example:
//
// itlib::static_vector<int, 4> myvec; // a static vector of size 0 and capacity 4
// myvec.resize(2); // vector is {0,0}
// myvec[1] = 11; // vector is {0,11}
// myvec.push_back(7); // vector is {0,11,7}
// myvec.insert(myvec.begin() + 1, 3); // vector is {0,3,11,7}
// myvec.push_back(5); // runtime error
// myvec.erase(myvec.begin());  // vector is {3,11,7}
// myvec.resize(5); // runtime error
//
//                  Configuration
//
// The library has two configuration options. They can be set as #define-s
// before including the header file, but it is recommended to change the code
// of the library itself with the values you want, especially if you include
// the library in many compilation units (as opposed to, say, a precompiled
// header or a central header).
//
//                  Config out of range error handling
//
// An out of range error is a runtime error which is triggered when the vector
// needs to be resized with a size higher than its capacity.
// For example: itlib::static_vector<int, 100> v(101);
//
// This is set by defining ITLIB_STATIC_VECTOR_ERROR_HANDLING to one of the
// following values:
// * ITLIB_STATIC_VECTOR_ERROR_HANDLING_NONE - no error handling. Crashes WILL
//      ensue if the error is triggered.
// * ITLIB_STATIC_VECTOR_ERROR_HANDLING_THROW - std::out_of_range is thrown.
// * ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT - asserions are triggered.
// * ITLIB_STATIC_VECTOR_ERROR_HANDLING_RESCUE - the error is ignored but
//      sanity is (somewhat) preserved. Functions which trigger the error will
//      just bail without changing the vector.
// * ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT_AND_THROW - combines assert and
//      throw to catch errors more easily in debug mode
// * ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT_AND_RESCUE - combines assert and
//      rescue to catch errors in debug mode and silently bail in release mode.
//
// To set this setting by editing the file change the line:
// ```
// #   define ITLIB_STATIC_VECTOR_ERROR_HANDLING ITLIB_STATIC_VECTOR_ERROR_HANDLING_THROW
// ```
// to the default setting of your choice
//
//                  Config bounds checks:
//
// By default bounds checks are made in debug mode (via an asser) when accessing
// elements (with `at` or `[]`). Iterators are not checked (yet...)
//
// To disable them, you can define ITLIB_STATIC_VECTOR_NO_DEBUG_BOUNDS_CHECK
// before including the header.
//
//
//                  MISC
//
// * There is an unused (and unusable) allocator class defined inside
// static_vector. It's point is to be a sham for templates which refer to
// container::allocator. It also allows it to work with itlib::flat_map
//
//
//                  TESTS
//
// You can find unit tests for static_vector in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <type_traits>
#include <cstddef>
#include <iterator>

#define ITLIB_STATIC_VECTOR_ERROR_HANDLING_NONE  0
#define ITLIB_STATIC_VECTOR_ERROR_HANDLING_THROW 1
#define ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT 2
#define ITLIB_STATIC_VECTOR_ERROR_HANDLING_RESCUE 3
#define ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT_AND_THROW 4
#define ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT_AND_RESCUE 5

#if !defined(ITLIB_STATIC_VECTOR_ERROR_HANDLING)
#   define ITLIB_STATIC_VECTOR_ERROR_HANDLING ITLIB_STATIC_VECTOR_ERROR_HANDLING_THROW
#endif

#if ITLIB_STATIC_VECTOR_ERROR_HANDLING == ITLIB_STATIC_VECTOR_ERROR_HANDLING_NONE
#   define I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(cond, rescue_return)
#elif ITLIB_STATIC_VECTOR_ERROR_HANDLING == ITLIB_STATIC_VECTOR_ERROR_HANDLING_THROW
#   include <stdexcept>
#   define I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(cond, rescue_return) if (cond) throw std::out_of_range("itlib::static_vector out of range")
#elif ITLIB_STATIC_VECTOR_ERROR_HANDLING == ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT
#   include <cassert>
#   define I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(cond, rescue_return) assert(!(cond) && "itlib::static_vector out of range")
#elif ITLIB_STATIC_VECTOR_ERROR_HANDLING == ITLIB_STATIC_VECTOR_ERROR_HANDLING_RESCUE
#   define I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(cond, rescue_return) if (cond) return rescue_return
#elif ITLIB_STATIC_VECTOR_ERROR_HANDLING == ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT_AND_THROW
#   include <stdexcept>
#   include <cassert>
#   define I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(cond, rescue_return) \
    do { if (cond) { assert(false && "itlib::static_vector out of range"); throw std::out_of_range("itlib::static_vector out of range"); } } while(false)
#elif ITLIB_STATIC_VECTOR_ERROR_HANDLING == ITLIB_STATIC_VECTOR_ERROR_HANDLING_ASSERT_AND_RESCUE
#   include <cassert>
#   define I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(cond, rescue_return) \
    do { if (cond) { assert(false && "itlib::static_vector out of range"); return rescue_return; } } while(false)
#else
#error "Unknown ITLIB_STATIC_VECTOR_ERRROR_HANDLING"
#endif


#if defined(ITLIB_STATIC_VECTOR_NO_DEBUG_BOUNDS_CHECK)
#   define I_ITLIB_STATIC_VECTOR_BOUNDS_CHECK(i)
#else
#   include <cassert>
#   define I_ITLIB_STATIC_VECTOR_BOUNDS_CHECK(i) assert((i) < this->size())
#endif

namespace itlib
{

template<typename T, size_t Capacity>
struct static_vector
{
    struct fake_allocator
    {
        typedef T value_type;
        typedef T* pointer;
        typedef const T* const_pointer;
    };

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
    using allocator_type = fake_allocator;

    static_vector() = default;

    explicit static_vector(size_t count)
    {
        resize(count);
    }

    static_vector(size_t count, const T& value)
    {
        assign_impl_val(count, value);
    }

    template <class InputIterator, typename = decltype(*std::declval<InputIterator>())>
    static_vector(InputIterator first, InputIterator last)
    {
        assign_impl_iter(first, last);
    }

    static_vector(std::initializer_list<T> l)
    {
        assign_impl_ilist(l);
    }

    static_vector(const static_vector& v)
    {
        for (const auto& i : v)
        {
            push_back(i);
        }
    }

    template <size_t CapacityB>
    static_vector(const static_vector<T, CapacityB>& v)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(v.size() > Capacity, );

        for (const auto& i : v)
        {
            push_back(i);
        }
    }

    static_vector(static_vector&& v) noexcept(std::is_nothrow_move_constructible<T>::value)
    {
        for (auto i = v.begin(); i != v.end(); ++i)
        {
            emplace_back(std::move(*i));
        }
        v.clear();
    }

    ~static_vector()
    {
        clear();
    }

    static_vector& operator=(const static_vector& v)
    {
        if (this == &v)
        {
            // prevent self usurp
            return *this;
        }

        clear();
        for (auto& elem : v)
        {
            push_back(elem);
        }

        return *this;
    }

    static_vector& operator=(static_vector&& v) noexcept(std::is_nothrow_move_assignable<T>::value)
    {
        clear();
        for (auto i = v.begin(); i!=v.end(); ++i)
        {
            emplace_back(std::move(*i));
        }

        v.clear();
        return *this;
    }

    void assign(size_type count, const T& value)
    {
        clear();
        assign_impl_val(count, value);
    }

    template <class InputIterator, typename = decltype(*std::declval<InputIterator>())>
    void assign(InputIterator first, InputIterator last)
    {
        clear();
        assign_impl_iter(first, last);
    }

    void assign(std::initializer_list<T> ilist)
    {
        clear();
        assign_impl_ilist(ilist);
    }

    const_reference at(size_type i) const
    {
        I_ITLIB_STATIC_VECTOR_BOUNDS_CHECK(i);
        return *reinterpret_cast<const T*>(m_data + i);
    }

    reference at(size_type i)
    {
        I_ITLIB_STATIC_VECTOR_BOUNDS_CHECK(i);
        return *reinterpret_cast<T*>(m_data + i);
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
        return reinterpret_cast<const T*>(m_data);
    }

    pointer data() noexcept
    {
        return reinterpret_cast<T*>(m_data);
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

    static constexpr size_t max_size() noexcept
    {
        return Capacity;
    }

    static constexpr size_t capacity() noexcept
    {
        return Capacity;
    }

    // modifiers
    void pop_back()
    {
        reinterpret_cast<const T*>(m_data + m_size - 1)->~T();
        --m_size;
    }

    void clear() noexcept
    {
        while (!empty())
        {
            pop_back();
        }
    }

    void push_back(const_reference v)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(size() >= Capacity, );

        ::new(m_data + m_size) T(v);
        ++m_size;
    }

    void push_back(T&& v)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(size() >= Capacity, );

        ::new(m_data + m_size) T(std::move(v));
        ++m_size;
    }

    template<typename... Args>
    reference emplace_back(Args&&... args)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(size() >= Capacity, );

        ::new(m_data + m_size) T(std::forward<Args>(args)...);
        ++m_size;
        return back();
    }

    iterator insert(iterator position, const value_type& val)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(size() >= Capacity, position);

        if (position == end())
        {
            emplace_back(val);
            return position;
        }

        emplace_back(std::move(back()));

        for (auto i = end() - 2; i != position; --i)
        {
            *i = std::move(*(i - 1));
        }

        *position = val;

        return position;
    }

    template<typename... Args>
    iterator emplace(iterator position, Args&&... args)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(size() >= Capacity, position);

        if (position == end())
        {
            emplace_back(std::forward<Args>(args)...);
            return position;
        }

        emplace_back(std::move(back()));

        for (auto i = end() - 2; i != position; --i)
        {
            *i = std::move(*(i - 1));
        }

        position->~T();
        ::new (position) T(std::forward<Args>(args)...);

        return position;
    }

    iterator erase(const_iterator position)
    {
        auto dist = position - begin();
        position->~T();

        for (auto i = begin() + dist + 1; i != end(); ++i)
        {
            *(i - 1) = std::move(*i);
        }

        resize(size() - 1);
        return begin() + dist;
    }

    void resize(size_type n)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(n > Capacity, );

        while (m_size > n)
        {
            pop_back();
        }

        while (n > m_size)
        {
            emplace_back();
        }
    }

    void swap(static_vector& v)
    {
        static_vector* longer;
        static_vector* shorter;

        if (v.m_size > m_size)
        {
            longer = &v;
            shorter = this;
        }
        else
        {
            longer = this;
            shorter = &v;
        }

        for (size_t i = 0; i < shorter->size(); ++i)
        {
            std::swap(shorter->at(i), longer->at(i));
        }

        auto short_size = shorter->m_size;

        for (size_t i = shorter->size(); i < longer->size(); ++i)
        {
            shorter->emplace_back(std::move(longer->at(i)));
            longer->at(i).~T();
        }

        longer->m_size = short_size;
    }

private:
    void assign_impl_val(size_t count, const T& value)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(count > Capacity, );

        for (size_t i = 0; i < count; ++i)
        {
            push_back(value);
        }
    }

    template <class InputIterator>
    void assign_impl_iter(InputIterator first, InputIterator last)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(long(last - first) > long(Capacity), );

        for (auto i = first; i != last; ++i)
        {
            push_back(*i);
        }
    }

    void assign_impl_ilist(std::initializer_list<T> l)
    {
        I_ITLIB_STATIC_VECTOR_OUT_OF_RANGE_IF(l.size() > Capacity, );

        for (auto&& i : l)
        {
            push_back(i);
        }
    }

    typename std::aligned_storage<sizeof(T), std::alignment_of<T>::value>::type m_data[Capacity];
    size_t m_size = 0;
};

template <typename T, size_t CapacityA, size_t CapacityB>
bool operator==(const static_vector<T, CapacityA>& a, const static_vector<T, CapacityB>& b)
{
    if (a.size() != b.size())
    {
        return false;
    }

    for (size_t i = 0; i < a.size(); ++i)
    {
        if (!(a[i] == b[i]))
            return false;
    }

    return true;
}

template <typename T, size_t CapacityA, size_t CapacityB>
bool operator!=(const static_vector<T, CapacityA>& a, const static_vector<T, CapacityB>& b)
{
    return !operator==(a, b);
}

}
