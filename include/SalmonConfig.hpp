/**
>HEADER
    Copyright (c) 2014, 2015, 2016, 2017 Rob Patro robp@cs.stonybrook.edu

    This file is part of Salmon.

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

#ifndef SALMON_CONFIG_HPP
#define SALMON_CONFIG_HPP

#include <string>

namespace salmon {
constexpr char majorVersion[] = "0";
constexpr char minorVersion[] = "14";
constexpr char patchVersion[] = "0";
constexpr char version[] = "0.14.0";
constexpr uint32_t indexVersion = 4;
constexpr char requiredQuasiIndexVersion[] = "q6";
} // namespace salmon

#endif // SALMON_CONFIG_HPP
