#include <stdint.h>

#include <algorithm>
#include <functional>

#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <jellyfish/misc.hpp>

namespace {
using testing::Types;

template <typename T>
class FloorLog2Test : public ::testing::Test {
public:
  T x;

  FloorLog2Test() : x(1) {}
};

typedef Types<uint32_t, uint64_t> Implementations;
TYPED_TEST_CASE(FloorLog2Test, Implementations);

TYPED_TEST(FloorLog2Test, FloorLog2) {
  unsigned int i = 0;

  for(i = 0; i < 8 * sizeof(this->x); i++) {
    ASSERT_EQ(i, jellyfish::floorLog2(this->x << i));
    ASSERT_EQ(i, jellyfish::floorLog2((this->x << i) | ((this->x << i) - 1)));
  }
}

using jellyfish::bitsize;
TEST(BitSize, Small) {
  for(int i = 1; i < 4098; ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    EXPECT_GE((1 << bitsize(i)) - 1, i);
    EXPECT_LE(1 << (bitsize(i) - 1), i);
  }
}

template<typename T>
int generic_popcount(T x) {
  static int counts[16] = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
  int res = 0;
  for( ; x; x >>= 4)
    res += counts[x & 0xf];
  return res;
}

TEST(Random, Bits) {
  uint64_t m        = 0;
  uint64_t m2       = 0;
  int      not_zero = 0;
  int      bits     = 0;
  for(int i = 0; i < 1024; ++i) {
    m        += jellyfish::random_bits(41);
    m2       += random_bits(41);
    not_zero += random_bits() > 0;
    bits     += generic_popcount(m);
  }
  std::cout << bits << "\n";
  EXPECT_LT((uint64_t)1 << 49, m); // Should be false with very low probability
  EXPECT_LT((uint64_t)1 << 49, m2); // Should be false with very low probability
  EXPECT_LT((int)(20.5 * 1024) / 2, bits);
  EXPECT_LT(512, not_zero);
}

TEST(BinarySearchFirst, Int) {
  static int size = 1024;

  for(int i = 0; i < size; ++i)
    EXPECT_EQ(i, *jellyfish::binary_search_first_false(jellyfish::pointer_integer<int>(0), jellyfish::pointer_integer<int>(size),
                                                       std::bind2nd(std::less<int>(), i)));
}

TEST(Slices, NonOverlapAll) {
  for(int iteration = 0; iteration < 100; ++iteration) {
    unsigned int nb_slices = random_bits(4) + 1;
    unsigned int size      = std::max(nb_slices, (unsigned int)random_bits(20));
    SCOPED_TRACE(::testing::Message() << "iteration:" << iteration
                 << " size:" << size << " nb_slices:" << nb_slices);

    unsigned int total = 0;
    unsigned int prev  = 0;
    for(unsigned int i = 0; i < nb_slices; ++i) {
      SCOPED_TRACE(::testing::Message() << "i:" << i);
      std::pair<unsigned int, unsigned int> b = jellyfish::slice(i, nb_slices, size);
      ASSERT_EQ(prev, b.first);
      ASSERT_GT(b.second, b.first);
      total += b.second - b.first;
      prev   = b.second;
    }
    ASSERT_EQ(size, total);
  }
}

TEST(QuoteArg, Nothing) {
  const std::string arg = "hello_world.1234_coucou/voila-bouh";
  EXPECT_EQ(arg, jellyfish::quote_arg(arg));
}
TEST(QuoteArg, Quote) {
  const std::string arg = "coucou voila";
  EXPECT_EQ("'" + arg + "'", jellyfish::quote_arg(arg));
}
TEST(QuoteArg, QuoteSlash) {
  const std::string arg = "coucou_'voila";
  EXPECT_EQ("'coucou_'\\''voila'", jellyfish::quote_arg(arg));
}
}
