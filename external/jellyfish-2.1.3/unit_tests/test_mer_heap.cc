#include <gtest/gtest.h>

#include <jellyfish/mer_dna.hpp>
#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/mer_heap.hpp>

namespace {
using jellyfish::mer_dna;

typedef jellyfish::large_hash::array<mer_dna>               large_array;
typedef large_array::region_iterator                        region_iterator;
typedef jellyfish::mer_heap::heap<mer_dna, region_iterator> mer_heap;

static const size_t   ary_size = 10000;
static const uint16_t mer_len  = 50;
static const size_t   nb_mers  = ary_size / 2;

class MerHeapTest : public ::testing::TestWithParam<size_t> {
protected:
  static void SetUpTestCase() {
    mer_dna::k(mer_len);
    shared_ary = new large_array(ary_size, mer_len * 2, 0, 63);
    mer_dna m;

    for(size_t i = 0; i < nb_mers; ++i) {
      m.randomize();
      bool   is_new;
      size_t id;
      shared_ary->set(m, &is_new, &id);
      EXPECT_TRUE(is_new); // Very small probability to fail
    }
  }

  static void TearDownTestCase() {
    delete shared_ary;
  }

  static large_array* shared_ary;
};
large_array* MerHeapTest::shared_ary = 0;

TEST_P(MerHeapTest, Order) {
  uint64_t           hash  = 0;
  int                count = 0;
  const large_array& ary   = *shared_ary;
  mer_dna            m;
  m.polyA();

  for(size_t slice = 0; slice < GetParam(); ++slice) {
    region_iterator rit = ary.iterator_slice<region_iterator>(slice, GetParam());
    mer_heap heap(ary.max_reprobe_offset());

    heap.fill(rit);

    while(!heap.is_empty()) {
      uint64_t nhash = ary.matrix().times(heap.head()->key_);
      EXPECT_LE(hash, nhash);
      if(nhash == hash)
        EXPECT_LT(m, heap.head()->key_);
      hash = nhash;
      m    = heap.head()->key_;

      heap.pop();
      ++count;

      if(rit.next())
        heap.push(rit);
    }
  }

  EXPECT_EQ(nb_mers, count);
}

INSTANTIATE_TEST_CASE_P(MerHeap, MerHeapTest, ::testing::Range((size_t)1, (size_t)11));

} // namespace {
