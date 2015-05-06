#include <gtest/gtest.h>
#include <stdlib.h>
#include <string.h>
#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <jellyfish/offsets_key_value.hpp>

namespace {
#ifndef UINT64_C
#define UINT64_C(x) ((uint64_t)x)
#endif

using namespace jellyfish;

// Tuple is {key_len, val_len, reprobe_len}. Actual reprobe value is computed from the reprobe_len.
class ComputeOffsetsTest : public ::testing::TestWithParam< ::std::tr1::tuple<int,int, int> >
{
public:
  Offsets<uint64_t> offsets;

  ComputeOffsetsTest() :
    offsets(::std::tr1::get<0>(GetParam()), ::std::tr1::get<1>(GetParam()), (1 << ::std::tr1::get<2>(GetParam())) - 2)
  { }

  ~ComputeOffsetsTest() { }
};

void test_key_offsets(const Offsets<uint64_t>::offset_t* it, unsigned int k_len, const char* message) {
  SCOPED_TRACE(message);

  unsigned int val_boff = (it->key.boff + k_len + 1 + (k_len - (bsizeof(uint64_t) - it->key.boff - 1)) / (bsizeof(uint64_t) - 1)) % bsizeof(uint64_t);
  if(it->key.sb_mask1) { // key between two words
    EXPECT_LT(sizeof(uint64_t), it->key.boff + k_len) <<
      ": Key not between words, should not have mask2";
    EXPECT_EQ(val_boff + (val_boff > 0), it->val.boff) <<
      ": Invalid value offset";
    if(it->val.boff == 0)
      EXPECT_EQ((uint64_t)0, it->key.sb_mask2);
  } else {
    if(bsizeof(uint64_t) >= it->key.boff + k_len) {
      EXPECT_EQ((it->key.boff + k_len) % bsizeof(uint64_t), it->val.boff) <<
        ": Invalid value offset";
    } else {
      EXPECT_TRUE(val_boff % bsizeof(uint64_t) == 0) <<
        ": Key between words, should have mask2";
    }
  }
  EXPECT_EQ((uint64_t)1 << (it->key.boff - 1), it->key.lb_mask) <<
    ": Invalid lb_mask";
  if(it->key.sb_mask1) {
    EXPECT_EQ(bsizeof(uint64_t) - it->key.boff - 1, it->key.shift) <<
      ": invalid key shift";
    EXPECT_EQ(it->key.boff == 1 ? 0 : (uint64_t)1 << (bsizeof(uint64_t) - it->key.boff + 1), 1 + (it->key.mask1 >> (it->key.boff - 1))) <<
      ": invalid key mask";
    EXPECT_EQ((uint64_t)1 << it->val.boff, 1 + it->key.mask2) <<
      ": invalid key mask2";
    EXPECT_EQ(bsizeof(uint64_t) - it->key.boff - 1, it->key.shift) <<
      ": invalid key shift";
    EXPECT_EQ((unsigned int)0, (k_len - (it->key.shift + it->key.cshift)) % (bsizeof(uint64_t) - 1)) <<
      ": invalid sum of shift and cshift";
    if(it->key.sb_mask2)
      EXPECT_EQ(1 + it->key.mask2, (uint64_t)1 << (it->key.cshift + 1)) <<
        ": invalid key cshift";
    else
      EXPECT_EQ((unsigned int)0, it->key.cshift);
    EXPECT_EQ((uint64_t)1 << 63, it->key.sb_mask1) <<
      ": invalid key key sb_mask1";
    EXPECT_EQ(k_len - (bsizeof(uint64_t) - 1 - it->key.boff) >= bsizeof(uint64_t) - 1, it->key.full_words);
  } else {
    EXPECT_EQ(0u, it->key.shift) <<
      ": key shift should be zero, no sb_mask";
    EXPECT_EQ((uint64_t)1 << (k_len + 1), 1 + (it->key.mask1 >> (it->key.boff - 1))) <<
      ": invalid key mask, no sb_mask";
    EXPECT_EQ(0u, it->key.cshift) <<
      ": key cshift should be zero, no sb_mask";
    EXPECT_EQ(0u, it->key.sb_mask1 | it->key.sb_mask2) <<
      ": set bit masks should be 0";
  }
}

void test_val_offsets(const Offsets<uint64_t>::offset_t* it, unsigned int v_len, const char* message) {
  SCOPED_TRACE(message);

  if(it->val.mask2) {
    EXPECT_TRUE(bsizeof(uint64_t) <= it->val.boff + v_len && (it->val.boff + v_len) % bsizeof(uint64_t) != 0) <<
      ": val not between words, should not have mask2";
    EXPECT_EQ(bsizeof(uint64_t) - it->val.boff, it->val.shift) <<
      ": invalid val shift";
    EXPECT_EQ((uint64_t)1 << (bsizeof(uint64_t) - it->val.boff), 1 + (it->val.mask1 >> it->val.boff)) <<
      ": invalid val mask1";
    EXPECT_EQ((uint64_t)1 << (it->val.boff + v_len - bsizeof(uint64_t)), 1 + it->val.mask2) <<
      ": invalid val mask2";
  } else {
    EXPECT_TRUE(bsizeof(uint64_t) > it->val.boff + v_len ||
                (bsizeof(uint64_t) <= it->val.boff + v_len && (it->val.boff + v_len) % bsizeof(uint64_t) == 0)) <<
      ": val between words, should have mask2";
    EXPECT_EQ(v_len, it->val.shift) <<
      ": invalid val shift";
    EXPECT_EQ(v_len == bsizeof(uint64_t) ? (uint64_t)0 : (uint64_t)1 << v_len, 1 + (it->val.mask1 >> it->val.boff)) <<
      ": invalid val mask1";
  }
  EXPECT_EQ(v_len, it->val.shift + it->val.cshift);
}

TEST_P(ComputeOffsetsTest, CheckCoherency) {
  const Offsets<uint64_t>::offset_t *it     = NULL, *pit = NULL;
  const Offsets<uint64_t>::offset_t *lit    = NULL, *lpit = NULL;
  unsigned int                             k_len  = ::std::tr1::get<0>(GetParam());
  unsigned int                             v_len  = ::std::tr1::get<1>(GetParam());
  unsigned int                             kv_len = k_len + v_len;
  unsigned int                             lk_len = ::std::tr1::get<2>(GetParam());
  unsigned int                             lv_len = std::min(kv_len - lk_len, (unsigned int)bsizeof(uint64_t));
  unsigned int                             i      = 0;

  EXPECT_EQ(lk_len, this->offsets.reprobe_len());
  EXPECT_EQ((uint64_t)1 << lk_len, this->offsets.reprobe_mask() + 1);
  EXPECT_EQ(lv_len, this->offsets.lval_len());


  for(i = 0; i < bsizeof(uint64_t); i++) {
    SCOPED_TRACE(::testing::Message() << "key:" << k_len << " val:" << v_len << " i:" << i 
                 << " lkey:" << lk_len << " lval:" << lv_len);
    this->offsets.word_offset(i, &it, &lit, NULL);
    if(i > 0) {
      this->offsets.word_offset(i-1, &pit, &lpit, NULL);

      EXPECT_GE(it->key.woff, pit->key.woff) << "Current key offset larger than previous";
      EXPECT_GE(it->key.woff, pit->val.woff) << "Current key offset larger than previous val offset";
      EXPECT_TRUE(it->key.woff != pit->key.woff || it->key.boff > pit->key.boff)
        << "Key bit offset did not increase within a word";
      EXPECT_EQ((pit->val.boff + v_len) % 64 + 1, it->key.boff) <<
        ": Invalid jump of bit offset";
    } // if(i > 0)
    EXPECT_GE(it->val.woff, it->key.woff) << "Val offset larger than key offset";
    EXPECT_TRUE(it->key.woff != it->val.woff || it->val.boff > it->key.boff)
      << "Val bit offset did not increase within a word";

    EXPECT_EQ(it->key.woff, lit->key.woff);
    EXPECT_EQ(it->key.boff, lit->key.boff);
    EXPECT_EQ(it->key.lb_mask, lit->key.lb_mask);

    test_key_offsets(it, k_len, "Key offsets");
    test_key_offsets(lit, lk_len, "Large key offsets");

    test_val_offsets(it, v_len, "Val offsets");
    test_val_offsets(lit, lv_len, "Large val offsets");

    size_t next_bit = (it->val.boff + v_len) % bsizeof(uint64_t);
    if(next_bit >= 62 || next_bit == 0)
      break;
  }

  EXPECT_EQ(i + 1, this->offsets.block_len());
}

INSTANTIATE_TEST_CASE_P(OffsetsTest, ComputeOffsetsTest, ::testing::Combine(::testing::Range(8, 192, 2), // Key lengths
                                                                            ::testing::Range(0, 64),    // Val lengths
                                                                            ::testing::Range(4, 9)      // Reprobe lengths
                                                                            ));
} // namespace {
