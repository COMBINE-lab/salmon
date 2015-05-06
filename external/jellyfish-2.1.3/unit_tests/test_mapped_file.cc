#include <iostream>
#include <fstream>
#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <jellyfish/mapped_file.hpp>

namespace {
using jellyfish::mapped_file;
TEST(MappedFile, CreateMove) {
  const char* mpt = "mapped_file_test";
  file_unlink file(mpt);
  std::string text = "Hello\nThis is a test\n";

  {
    std::ofstream fd(mpt);
    fd << text;
  }

  mapped_file mf(mpt);
  EXPECT_STREQ(mpt, mf.path().c_str());
  ASSERT_NE((char*)0, mf.base());
  ASSERT_EQ(text.size(), mf.length());
  EXPECT_EQ(text, std::string(mf.base(), mf.length()));

  mapped_file mf2 = std::move(mf);
  EXPECT_STREQ(mpt, mf2.path().c_str());
  ASSERT_NE((char*)0, mf2.base());
  ASSERT_EQ(text.size(), mf2.length());
  EXPECT_EQ(text, std::string(mf2.base(), mf2.length()));
  ASSERT_EQ((char*)0, mf.base());

  mapped_file mf3(mpt);
  mf2 = std::move(mf3);
  EXPECT_STREQ(mpt, mf2.path().c_str());
  ASSERT_NE((char*)0, mf2.base());
  ASSERT_EQ(text.size(), mf2.length());
  EXPECT_EQ(text, std::string(mf2.base(), mf2.length()));
  ASSERT_EQ((char*)0, mf.base());
}

// Gtest and newer compilers seem to have a problem with EXPECT_THROW
#if !defined(__clang__) && (!defined(GTEST_GCC_VER_) || GTEST_GCC_VER_ < 40800)
TEST(MappedFile, Fail) {
  const char* bad_file = "/doesntexistsforsure/thatwouldbecrazy!";
  EXPECT_THROW(mapped_file mf(bad_file), jellyfish::mapped_file::ErrorMMap);

  mapped_file mf;
  EXPECT_THROW(mf.map(bad_file), jellyfish::mapped_file::ErrorMMap);
  EXPECT_EQ((char*)0, mf.base());
  EXPECT_THROW(mf.map(-1), jellyfish::mapped_file::ErrorMMap);
  EXPECT_EQ((char*)0, mf.base());
}
#endif
}
