// itlib-qalgorithm v1.01
//
// Wrappers of <algorithm> algorithms for entire containers
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
//  1.01 (2020-12-29) Added erase functions
//  1.00 (2020-12-28) First pulic release
//
//
//                  DOCUMENTATION
//
// Simply include this file wherever you need.
// It defines the following algorithms:
//
// * qfind - wraps std::find
// * qfind_if - wraps std::find_if
// * pfind - wraps std::find, returns a raw pointer to the element or nullptr if the element wasn't found
// * pfind_if - wraps std::find_if, returns a raw pointer to the element or nullptr if the element wasn't found
// * bool erase_first(container, value) - erase the first element equal to value. returns true if something was erased
// * bool erase_first_if(container, pred) - erase the first elemenf which matches pred. returns true if something was erased
// * size_t erase_all(container, value) - erases all elements equal to value, returns number of elements erased
// * size_t erase_all_if(container, value) - erases all elements which match pred, returns number of elements erased
//
//
//                  TESTS
//
// You can find unit tests for qalgorithm in its official repo:
// https://github.com/iboB/itlib/blob/master/test/
//
#pragma once

#include <utility>
#include <algorithm>
#include <type_traits>

namespace itlib
{

namespace impl
{
// get the appropriate return type: iterator for non-const containers and const_iterator for const containers
template <typename Container>
struct iterator_t
{
    using type = typename std::conditional<std::is_const<Container>::value, typename Container::const_iterator, typename Container::iterator>::type;
};
template <typename Container>
struct pointer_t
{
    using type = typename std::conditional<std::is_const<Container>::value, typename Container::value_type const*, typename Container::value_type*>::type;
};
}

template <typename Container, typename Value>
typename impl::iterator_t<Container>::type qfind(Container& c, const Value& val)
{
    return std::find(c.begin(), c.end(), val);
}

template <typename Container, typename Value>
typename impl::pointer_t<Container>::type pfind(Container& c, const Value& val)
{
    auto f = std::find(c.begin(), c.end(), val);
    if (f == c.end()) return nullptr;
    return &(*f);
}

template <typename Container, typename Pred>
typename impl::iterator_t<Container>::type qfind_if(Container& c, Pred&& pred)
{
    return std::find_if(c.begin(), c.end(), std::forward<Pred>(pred));
}

template <typename Container, typename Pred>
typename impl::pointer_t<Container>::type pfind_if(Container& c, Pred&& pred)
{
    auto f = std::find_if(c.begin(), c.end(), std::forward<Pred>(pred));
    if (f == c.end()) return nullptr;
    return &(*f);
}

template <typename Container, typename Value>
bool erase_first(Container& c, const Value& val)
{
    auto f = qfind(c, val);
    if (f == c.end()) return false;
    c.erase(f);
    return true;
}

template <typename Container, typename Pred>
bool erase_first_if(Container& c, Pred&& pred)
{
    auto f = qfind_if(c, std::forward<Pred>(pred));
    if (f == c.end()) return false;
    c.erase(f);
    return true;
}

template <typename Container, typename Value>
typename Container::size_type erase_all(Container& c, const Value& val)
{
    auto newend = std::remove(c.begin(), c.end(), val);
    auto ret = c.end() - newend;
    c.erase(newend, c.end());
    return ret;
}

template <typename Container, typename Pred>
typename Container::size_type erase_all_if(Container& c, Pred&& pred)
{
    auto newend = std::remove_if(c.begin(), c.end(), std::forward<Pred>(pred));
    auto ret = c.end() - newend;
    c.erase(newend, c.end());
    return ret;
}

}
