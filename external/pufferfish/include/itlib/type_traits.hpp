// itlib-type_traits v1.00
//
// Additional helper type traits extending the standard <type_traits>
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
//  1.00 (2020-12-28) First pulic release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the following classes:
//
// * is_instantiation_of<Template, Type> - checks whether a type is an
//   instantiation of a template
//      Example:
//      using MyVec = std::vector<int>;
//      static_assert(is_instantiation_of<std::vector, MyVec>::value,
//          "MyVec must be a std::vector");
//
// With C++17 all traits have a _v template constant
//
//
//                  TESTS
//
// You can find unit tests for type_traits in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <type_traits>

namespace itlib
{

template <template <typename...> class, typename...>
struct is_instantiation_of : public std::false_type {};

template <template <typename...> class Template, typename... TArgs>
struct is_instantiation_of<Template, Template<TArgs...>> : public std::true_type {};

#if __cplusplus >= 201700
template <template <typename...> class Template, typename Type>
inline constexpr bool is_instantiation_of_v = is_instantiation_of<Template, Type>::value;
#endif

}
