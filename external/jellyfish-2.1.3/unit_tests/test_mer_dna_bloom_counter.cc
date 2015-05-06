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
#include <algorithm>
#include <utility>
#include <set>
#include <string>
#include <fstream>
#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>

namespace {
using jellyfish::mer_dna;

static const size_t nb_inserts = 10000;
static const double error_rate = 0.001;

template<typename T>
class MerDnaBloomTest : public ::testing::Test { };

struct TestBloomCounter {
  typedef jellyfish::mer_dna_bloom_counter bloom_type;
  typedef jellyfish::mer_dna_bloom_counter_file file_type;
  static const unsigned int threshold_twice = 2; // Bloom counter counts up to 2.
};
struct TestBloomFilter {
  typedef jellyfish::mer_dna_bloom_filter bloom_type;
  typedef jellyfish::mer_dna_bloom_filter_file file_type;
  static const unsigned int threshold_twice = 1; // Bloom filter counts up to 1.
};

typedef ::testing::Types<TestBloomCounter, TestBloomFilter> TestBloomCounterTypes;
TYPED_TEST_CASE(MerDnaBloomTest, TestBloomCounterTypes);


TYPED_TEST(MerDnaBloomTest, FalsePositive) {
  mer_dna::k(50);
  std::set<mer_dna> mer_set;
  typename TypeParam::bloom_type bc(error_rate, nb_inserts);

  size_t collisions2 = 0;
  size_t collisions3 = 0;

  // Insert once nb_inserts. Insert twice the first half
  {
    // First insertion
    size_t nb_collisions   = 0;
    mer_dna m;
    for(size_t i = 0; i < nb_inserts; ++i) {
      m.randomize();
      mer_set.insert(m);
      nb_collisions += bc.insert(m) > 0;
    }
    EXPECT_GT(error_rate * nb_inserts, nb_collisions);
  }

  // Second insertion
  {
    size_t nb_collisions = 0;
    size_t nb_errors     = 0;
    auto it = mer_set.cbegin();
    for(size_t i =0; i < nb_inserts / 2; ++i, ++it) {
      unsigned int oc = bc.insert(*it);
      nb_collisions += oc > 1;
      nb_errors += oc < 1;
    }
    EXPECT_GT(2 * error_rate * nb_inserts, nb_collisions);
    EXPECT_EQ((size_t)0, nb_errors);
  }

  // Write to file and reload two different ways
  file_unlink f("bloom_file");
  {
    std::ofstream out(f.path.c_str());
    bc.write_bits(out);
    EXPECT_TRUE(out.good());
    EXPECT_EQ(bc.nb_bytes(), out.tellp());
  }
  std::ifstream in(f.path.c_str());
  typename TypeParam::bloom_type bc_read(bc.m(), bc.k(), in, bc.hash_functions());
  EXPECT_EQ(bc.nb_bytes(), in.tellg());
  in.close();
  typename TypeParam::file_type bc_map(bc.m(), bc.k(), f.path.c_str(), bc.hash_functions());
  EXPECT_EQ(bc.m(), bc_read.m());
  EXPECT_EQ(bc.k(), bc_read.k());
  EXPECT_EQ(bc.m(), bc_map.m());
  EXPECT_EQ(bc.k(), bc_map.k());

  // Check known mers
  {
    size_t nb_collisions = 0;
    size_t nb_errors     = 0;
    auto it = mer_set.cbegin();
    for(size_t i = 0; i < nb_inserts; ++i, ++it) {
      unsigned int check = bc.check(*it);
      EXPECT_EQ(check, bc_read.check(*it));
      EXPECT_EQ(check, bc_map.check(*it));
      if(i < nb_inserts / 2) {
        nb_errors += check < TypeParam::threshold_twice;
      } else {
        nb_errors += check < 1;
        nb_collisions += check > 1;
      }
    }
    EXPECT_EQ((size_t)0, nb_errors);
    EXPECT_GT(2 * error_rate * nb_inserts, nb_collisions);
  }

  // Check unknown mers
  {
    size_t nb_collisions = 0;
    mer_dna m;
    for(size_t i = 0; i < nb_inserts; ++i) {
      m.randomize();
      unsigned int check = bc.check(m);
      EXPECT_EQ(check, bc_read.check(m));
      EXPECT_EQ(check, bc_map.check(m));
      nb_collisions += check > 0;
    }
    EXPECT_GT(2 * error_rate * nb_inserts, nb_collisions);
  }
}

TYPED_TEST(MerDnaBloomTest, Move) {
  mer_dna::k(100);
  typename TypeParam::bloom_type bc(error_rate, nb_inserts);
  const unsigned long            k = bc.k();
  const size_t                   m = bc.m();
  typename TypeParam::bloom_type bc2(std::move(bc));
  EXPECT_EQ(k, bc2.k());
  EXPECT_EQ(m, bc.m());
}

}
