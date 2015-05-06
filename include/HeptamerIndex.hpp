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


#ifndef __HEPTAMER_INDEX_HPP__
#define __HEPTAMER_INDEX_HPP__

#include <atomic>
#include <vector>
#include <array>

class HeptamerIndex {
  using AtomicCount = std::atomic<uint64_t>;
public:
  explicit HeptamerIndex();
  std::size_t index(uint64_t heptamer);
private:
  constexpr static uint32_t PossibleHeptamers = 16384;
  const std::array<uint32_t, 7> mult_{{1,4,16,64,256,1024,4096}};
  std::vector<AtomicCount> heptamers_;
};

#endif // __HEPTAMER_INDEX_HPP__
