// itlib-expected v1.01
//
// A union-type of a value and an error
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2021 Borislav Stanimirov
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
//  1.01 (2021-09-27) Fixed value_or which could return a ref to temporary
//  1.00 (2021-09-26) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::expected, which is a union type of a value and
// an error.
//
// It is somewhat similar to std::optional. In a way you can think of
// std::optional<T> as an itlib::expected<T, bool> WITH THE NOTABLE DIFFERENCE
// that a default-constructed expected is truthy (it invokes the default
// contstructor) of T
//
// It is also similar to the the proposed std::expected
// http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2017/p0323r3.pdf
// however it is not the same. itlib::expected does 90% of the job of the
// proposed std::expected with 10% of the code.
//
// Notably itlib::expected is a non-copyable type. It is simply too much effort
// to take care of the case where a copy-constructor throws and copying it
// seems like a very, very rare need.
//
// In most cases one would return an itlib::expected from a function and the
// lifetime of the object will span only the caller's scope.
//
//                  Basics
//
// itlib::expected has two template arguments:
// * an "expected" or "value" type T: typedef-ed as value_type
// * an "unexpected" or "error" type E: typedef-ed as error_type
//
// It is a union type: It holds either a value or an error.
//
// It is a boolean type: It is truthy if it holds a value and falsy if it holds
// an error.
//
// It is an accessor type: It has dereference operators (* and ->) which lead
// to the value type inside. Invoking them if itlib::expected holds an error
// leads to undefined behavior.
//
// It is a non-copyable type: Only move ops are allowed, even if the value or
// error inside are copyable.
//
// Ways to construct it with an "expected" type:
// * The default constructor. It will in turn invoke the default constructor
// of value_type.
// * With an existing variable of value_type which can be copied or moved
// inside. itlib::expected has a implicit constructor from T&&, so you can
// simply return a T from a function which returns itlib::expected.
//
// Ways to costruct with an "unexpected" type
// * Use `unexpected()`. It can be use to initialize itlib::expected with the
// default constructor of error_type
// * Use `unexpected(E&&)`. It can construct an itlib::expected with an
// existing error type variable
//
//                  Example
//
//    enum class int_error { division_by_zero, signed_overflow };
//
//    itlib::expected<int, int_error> divide(int a, int b) {
//        if (b == 0) return unexpected(int_error::division_by_zero);
//        return a/b;
//    }
//    ...
//    auto res = divide(x, y);
//    if (!res) cerr << int_error_to_string(res.error()) << '\n';
//    else cout << *res << '\n';
//
//                  Reference
//
// Helpers
// * unexpected() - free function used to create an itlib::expected with a
//   default-constructed error_type
// * unexpected(E&&) - free function used to create an itlib::expected with a
//   error_type from E
// * unexpected_t<E> - type which can be used to create an itlib::expected with
//   any value_type and an error which can be constructed from E
// * unexpected_t<void> - type which can be used to create an itlib::expected
//   with any value_type and any error_type
// Main class: expected<T, E> - union type of a value T and error E
// * expected() - constructs with a default-constructed value T
// * expected(T&& t) - constructs with a value T from t
// * expected(unexpected_t<E>) - constructs with an error
// * expected(expected&&) - move ctor
// * operator=(expected&&) - move assignment
// Boolean interface:
// * has_value() - true if it holds a value
// * has_error() - true if it holds an error
// * operator bool() = has_value()
// Get value:
// * value() - returns value
// * value_or(T t) - returns value if truthy or t if falsy
// * operator* - returns value
// * operator-> - returns value
// Get error:
// * error() - return error
//
//                  TESTS
//
// You can find unit tests for expected in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <cassert>
#include <utility>
#include <new>

namespace itlib
{

template <typename E>
class unexpected_t
{
public:
    explicit unexpected_t(E&& e) : m_error(std::forward<E>(e)) {}

private:
    template<typename T, typename E2>
    friend class expected;

    E m_error;
};

template <typename E>
unexpected_t<E> unexpected(E&& e)
{
    return unexpected_t<E>(std::forward<E>(e));
}

template <>
class unexpected_t<void> {};

inline unexpected_t<void> unexpected() noexcept { return {}; }

template <typename T, typename E>
class expected
{
public:
    using value_type = T;
    using error_type = E;

    expected() : m_value(), m_has_value(true) {}
    expected(T&& t) : m_value(std::forward<T>(t)), m_has_value(true) {}

    template <typename E2>
    expected(unexpected_t<E2>&& u) : m_error(std::move(u.m_error)), m_has_value(false) {}

    expected(unexpected_t<void>) : m_error(), m_has_value(false) {}

    // do not copy
    expected(const expected&) = delete;
    expected& operator=(const expected&) = delete;

    // do move
    expected(expected&& other) noexcept
        : m_has_value(other.has_value())
    {
        if (m_has_value)
        {
            new (&m_value) T(std::move(other.m_value));
        }
        else
        {
            new (&m_error) E(std::move(other.m_error));
        }
    }

    expected& operator=(expected&& other) noexcept
    {
        if (m_has_value && other.has_value())
        {
            m_value = std::move(other.m_value);
        }
        else if(m_has_value && !other.has_value())
        {
            m_has_value = false;
            m_value.~T();
            ::new (&m_error) E(std::move(other.m_error));
        }
        else if(!m_has_value && other.has_value())
        {
            m_has_value = true;
            m_error.~E();
            ::new (&m_value) T(std::move(other.m_value));
        }
        else
        {
            m_error = std::move(other.m_error);
        }
        return *this;
    }

    ~expected()
    {
        if (m_has_value)
        {
            m_value.~T();
        }
        else
        {
            m_error.~E();
        }
    }

    // bool interface
    bool has_value() const { return m_has_value; }
    bool has_error() const { return !m_has_value; }
    explicit operator bool() const { return m_has_value; }

    // value getters
    T& value() &
    {
        assert(has_value());
        return m_value;
    }

    const T& value() const &
    {
        assert(has_value());
        return m_value;
    }

    T&& value() &&
    {
        assert(has_value());
        return std::move(m_value);
    }

    T& operator*() & { return value(); }
    const T& operator*() const & { return value(); }
    T&& operator*() && { return std::move(value()); }

    template <typename U>
    T value_or(U&& v) const & { return has_value() ? value() : std::forward<U>(v); }
    template <typename U>
    T value_or(U&& v) && { return has_value() ? std::move(value()) : std::forward<U>(v); }

    T* operator->() { return &value(); }
    const T* operator->() const { return &value(); }

    // error getters

    E& error() &
    {
        assert(has_error());
        return m_error;
    }

    const E& error() const &
    {
        assert(has_error());
        return m_error;
    }

    E&& error() &&
    {
        assert(has_error());
        return std::move(m_error);
    }

private:
    union
    {
        value_type m_value;
        error_type m_error;
    };
    bool m_has_value;
};

}
