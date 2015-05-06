/*  This file is part of Jellyfish.

    Jellyfish is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Jellyfish is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Jellyfish.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <unit_tests/test_main_cmdline.hpp>
#include <unit_tests/test_main.hpp>
#include <jellyfish/backtrace.hpp>

template<long int n>
struct floorLog2 {
  static const int val = floorLog2<n / 2>::val + 1;
};
template<>
struct floorLog2<1> {
  static const int val = 0;
};

// Return length random bits. 1 <= length <= 64
uint64_t random_bits(int length) {
  uint64_t res = 0;
  for(int i = 0; i < length; i += floorLog2<RAND_MAX>::val)
    res ^= (uint64_t)random() << i;
  return res & ((uint64_t)-1 >> (64 - length));
}


int main(int argc, char *argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  cmdline_parse args(argc, argv);

  unsigned int seed;
  if(args.seed_given) {
    seed = args.seed_arg;
  } else {
    std::ifstream urandom("/dev/urandom");
    urandom.read((char*)&seed, sizeof(seed));
    if(!urandom.good()) {
      std::cerr << "Failed to read random seed" << std::endl;
      return 1;
    }
  }
  if(args.backtrace_flag) {
    show_backtrace();
    setenv("GTEST_CATCH_EXCEPTIONS", "0", 1);
  }

  std::cout << "Using random seed " << seed << std::endl;
  srandom(seed);

  return RUN_ALL_TESTS();
}
