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

#ifndef __SCOPED_TIMER_HPP__
#define __SCOPED_TIMER_HPP__
// from https://gist.github.com/justgord/4482447
#include <chrono>
#include <iostream>

struct ScopedTimer {
  std::chrono::high_resolution_clock::time_point t0;

  ScopedTimer(bool print = true)
      : t0(std::chrono::high_resolution_clock::now()), print_(print) {}
  ~ScopedTimer(void) {
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedSec = t1 - t0;
    if (print_) {
      std::cerr << "Elapsed time: " << elapsedSec.count() << "s\n";
    }
  }

private:
  bool print_;
};

#endif //__SCOPED_TIMER_HPP__
