/**
>HEADER
    Copyright (c) 2014-2024 Rob Patro rob@cs.umd.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#ifndef SALMON_CONFIG_HPP
#define SALMON_CONFIG_HPP

#include <cstdint>
#include <string>

namespace salmon {
constexpr char majorVersion[] = "1";
constexpr char minorVersion[] = "10";
constexpr char patchVersion[] = "3";
constexpr char version[] = "1.10.3";
constexpr uint32_t indexVersion = 5;
constexpr char requiredQuasiIndexVersion[] = "p7";
} // namespace salmon

#endif // SALMON_CONFIG_HPP
