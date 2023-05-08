// itlib-sentry v1.01
//
// A sentry which invokes a function object when destroyed
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2020-2022 Borislav Stanimirov
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
//  1.01 (2022-01-13) [[nodicard]] of type if compiled with C++17
//  1.00 (2020-10-15) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the class itlib::sentry. The sentry invokes a provided function
// object (lambda, free function, object with operator()) when destroyed. It's
// useful when you want to invoke some code when exiting a scope, which has
// many exit points.
//
// To create a sentry with C++11, you can use itlib::make_sentry(func), which
// will properly construct the sentry.
//
// With C++17 and template argument deduction you can just create a local
// sentry object with the given function
//
// C++11 example:
//
// void foo(int n);
// {
//     itlib::sentry at_exit = itlib::make_sentry([]() { cout << "exiting foo"; });
//     if (n==0) return;
//     for (int i=0; i<n; ++i) bar();
// }
//
// C++17 example:
//
// void foo(int n);
// {
//     itlib::sentry at_exit([]() { cout << "exiting foo"; });
//     if (n==0) return;
//     for (int i=0; i<n; ++i) bar();
// }
//
//
//                  TESTS
//
// You can find unit tests for sentry in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <utility>

#if !defined(ITLIB_NODISCARD)
#   if __cplusplus >= 201700
#       define ITLIB_NODISCARD [[nodiscard]]
#   else
#       define ITLIB_NODISCARD
#   endif
#endif

namespace itlib
{

template <typename Func>
class ITLIB_NODISCARD sentry
{
public:
    explicit sentry(Func&& atexit) : m_func(std::forward<Func>(atexit)) {}

    sentry(const sentry&) = delete;
    sentry& operator=(const sentry&) = delete;

    sentry& operator=(sentry&&) = delete;

private:
    Func m_func;

public:
#if __cplusplus >= 201700
    // c++ 17 has guaranteed copy-elision, so we can afford to do this
    sentry(sentry&&) = delete;
    ~sentry() { m_func(); }
#else
    sentry(sentry&& other) noexcept
        : m_func(std::move(other.m_func))
        , m_has_func(other.m_has_func)
    {
        other.m_has_func = false;
    }
    ~sentry() { if (m_has_func) m_func(); }
private:
    bool m_has_func = true;
#endif
};

template <typename Func>
sentry<Func> make_sentry(Func&& func)
{
    return sentry<Func>(std::forward<Func>(func));
}

}
