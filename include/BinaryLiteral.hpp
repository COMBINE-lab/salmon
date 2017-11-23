/**
>HEADER
    Copyright (c) 2013 Rob Patro robp@cs.cmu.edu

    This file is part of Sailfish.

    Sailfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Sailfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Sailfish.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#ifndef __BINARY_LITERAL_HPP__
#define __BINARY_LITERAL_HPP__

//===============================================================

//
// binary_literal_impl represents a compile-time function that
// computes the unsigned long long int from a list of characters
// Digits that MUST be composed of '0' or '1'...
//
template <char... Digits> struct binary_literal_impl;

// If the next digit is zero, then compute the rest...
template <char... Digits> struct binary_literal_impl<'0', Digits...> {
  static constexpr unsigned long long to_ulonglong() {
    return binary_literal_impl<Digits...>::to_ulonglong();
  }
};

// If the next digit is one, then shift 1 and compute the rest...
template <char... Digits> struct binary_literal_impl<'1', Digits...> {
  static constexpr unsigned long long to_ulonglong() {
    return (1UL << sizeof...(Digits)) |
           binary_literal_impl<Digits...>::to_ulonglong();
  }
};

// Base case: No digits, so return 0...
template <> struct binary_literal_impl<> {
  static constexpr unsigned long long to_ulonglong() { return 0; }
};

//===============================================================

template <char... Digits> constexpr unsigned long long operator"" _binary() {
  return binary_literal_impl<Digits...>::to_ulonglong();
}

//===============================================================

#endif // __BINARY_LITERAL_HPP__
