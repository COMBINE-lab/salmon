//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __PAIRED_ALIGNMENT_FORMATTER_HPP__
#define __PAIRED_ALIGNMENT_FORMATTER_HPP__

//#include "RapMapUtils.hpp"
#include "Util.hpp"

template <typename IndexPtrT>
struct PairedAlignmentFormatter {
  PairedAlignmentFormatter(IndexPtrT indexIn) : index(indexIn),
                                                read1Temp(1000, 'A'),
                                                qual1Temp(1000, '~'),
                                                read2Temp(1000, 'A'),
                                                qual2Temp(1000, '~'),
                                                cigarStr1(buff1, 1000),
                                                cigarStr2(buff2, 1000),
                                                use_qualities(false) {
  }

  bool enable_qualities() {
    bool prev = use_qualities;
    use_qualities = true;
    return prev;
  }

  bool disable_qualities() {
    bool prev = use_qualities;
    use_qualities = false;
    return prev;
  }

  // Data members
  IndexPtrT index;
  std::string read1Temp;
  std::string qual1Temp;
  std::string read2Temp;
  std::string qual2Temp;
  char buff1[1000];
  char buff2[1000];
  pufferfish::util::FixedWriter cigarStr1;
  pufferfish::util::FixedWriter cigarStr2;
  bool use_qualities;
  std::string empty_qual{"*"};
};

#endif //__PAIR_ALIGNMENT_FORMATTER_HPP__
