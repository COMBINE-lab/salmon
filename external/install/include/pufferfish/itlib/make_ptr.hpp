// itlib-make-ptr v1.00
//
// Helper functions for making std::shared_ptr and std::unique_ptr
//
// SPDX-License-Identifier: MIT
// MIT License:
// Copyright(c) 2020 Borislav Stanimirov
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
//  1.00 (2020-10-15) Initial release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the following functions:
//
// * itlib::make_shared(T&& arg)
//   Create a std::shared_ptr<T> by invoking std::make_shared by forwarding arg
//   to T's constructor. Thus allowing to make shared pointers which hold a
//   a copy of arg, or in case arg is an rvalue, it will me moved to the pointer.
//   This allowing to to save a retype of the name in case you want to copy or
//   move a value into a new std::shared_ptr
//
// * itlib::make_unique(T&& arg)
//   The same as itlib::make_shared but for std::unique_ptr
//   It's also written for C++11, so you don't need to enable C++14 to include
//   this header. However it's not a full substitution for std::make_unique
//
// Example:
//
// my<complex, template, type> val;
//
// // not nice
// auto ptr = std::make_shared<my<complex, template, type>>(std::move(val));
//
// // nice
// auto p1 = itlib::make_shaerd(val); // copy val into p1
// auto p2 = itlib::make_shaerd(std::move(val)); // move val into p2
//
//
//                  TESTS
//
// You can find unit tests for make-ptr in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <memory>
#include <type_traits>

namespace itlib
{

template <typename T>
auto make_shared(T&& t) -> std::shared_ptr<typename std::remove_reference<T>::type>
{
    return std::make_shared<typename std::remove_reference<T>::type>(std::forward<T>(t));
}

template <typename T>
auto make_unique(T&& t) -> std::unique_ptr<typename std::remove_reference<T>::type>
{
    using RRT = typename std::remove_reference<T>::type;
    return  std::unique_ptr<RRT>(new RRT(std::forward<T>(t)));
}

}
