#include <gtest/gtest.h>
#include <jellyfish/allocators_mmap.hpp>

namespace {
TEST(AllocMmap, Simple) {
  static const size_t size = 4096;
  allocators::mmap mem(size);

  EXPECT_EQ(size, mem.get_size());
  EXPECT_NE((void*)0, mem.get_ptr());
  char* ptr = (char*)mem.get_ptr();
  // char c = ptr[size];
  // EXPECT_EQ((char)0, c);
}
}
