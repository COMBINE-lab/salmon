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

#ifndef COMMON_TYPES_HPP
#define COMMON_TYPES_HPP

#include <array>
#include <cstdint>

namespace Sailfish {
/*
struct Kmer {
    uint64_t kmer;
};

struct Index {
    uint64_t index;
};
*/
struct TranscriptFeatures {
  std::string name;
  size_t length;
  double gcContent;
  std::array<uint64_t, 16> diNucleotides;
};
} // namespace Sailfish

#endif // COMMON_TYPES_HPP
