#include <gtest/gtest.h>
#include <config.h>

#ifdef HAVE_INT128
#include <limits>
#include <jellyfish/int128.hpp>

namespace {
TEST(Int128, Specialized) {
  // This is supposed to be true, whether it is specialized by gcc or
  // in int128.hpp above.
  EXPECT_TRUE(std::numeric_limits<__int128>::is_specialized);
  EXPECT_TRUE(std::numeric_limits<unsigned __int128>::is_specialized);
}

TEST(Int128, Max) {
  static const __int128 sone = 1;
  static const unsigned __int128 uone = 1;
  __int128 sm = std::numeric_limits<__int128>::max();
  unsigned __int128 um = std::numeric_limits<unsigned __int128>::max();

  for(int i = 0; i < 127; ++i, sm >>= 1, um >>= 1) {
    EXPECT_EQ(1, (int)(sm & sone));
    EXPECT_EQ(1, (int)(um & uone));
    EXPECT_EQ(1, (int)((std::numeric_limits<__int128>::max() >> i) & sone));
    EXPECT_EQ(1, (int)((std::numeric_limits<unsigned __int128>::max() >> i) & uone));
  }

  EXPECT_EQ(0, (int)sm);
  EXPECT_EQ(1, (int)um);
}
} // namespace

#endif // HAVE_INT128
