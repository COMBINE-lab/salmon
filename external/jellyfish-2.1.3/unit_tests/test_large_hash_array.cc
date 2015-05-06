#include <sys/types.h>
#include <signal.h>
#include <unistd.h>

#include <map>
#include <vector>
#include <limits>

#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>

#include <jellyfish/large_hash_array.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/atomic_gcc.hpp>
#include <jellyfish/allocators_mmap.hpp>

void PrintTo(jellyfish::mer_dna& m, ::std::ostream* os) {
  *os << m.to_str();
}

namespace {
typedef jellyfish::large_hash::array<jellyfish::mer_dna> large_array;
typedef std::map<jellyfish::mer_dna, uint64_t> mer_map;
typedef std::set<jellyfish::mer_dna> mer_set;

using jellyfish::RectangularBinaryMatrix;
using jellyfish::mer_dna;
using std::numeric_limits;

typedef large_array::iterator stl_iterator;
typedef large_array::eager_iterator eager_iterator;
typedef large_array::lazy_iterator lazy_iterator;
typedef large_array::region_iterator region_iterator;

// Tuple is {key_len, val_len, reprobe_len}.
class HashArray : public ::testing::TestWithParam< ::std::tr1::tuple<int,int, int> >
{
public:
  static const size_t ary_lsize = 10;
  static const size_t ary_size = (size_t)1 << ary_lsize;
  static const size_t ary_size_mask = ary_size - 1;
  const int           key_len, val_len, reprobe_len, reprobe_limit;
  large_array         ary;

  HashArray() :
    key_len(::std::tr1::get<0>(GetParam())),
    val_len(::std::tr1::get<1>(GetParam())),
    reprobe_len(::std::tr1::get<2>(GetParam())),
    reprobe_limit((1 << reprobe_len) - 2),
    ary(ary_size, key_len, val_len, reprobe_limit)
  { }

  void SetUp() {
    jellyfish::mer_dna::k(key_len / 2);
  }

  ~HashArray() { }
};

TEST_P(HashArray, OneElement) {
  mer_dna m, m2, get_mer;

  SCOPED_TRACE(::testing::Message() << "key_len:" << key_len << " val_len:" << val_len << " reprobe:" << reprobe_limit);

  EXPECT_EQ((unsigned int)ary_lsize, ary.matrix().r());
  EXPECT_EQ((unsigned int)key_len, ary.matrix().c());

  size_t start_pos = random() % (ary_size - bsizeof(uint64_t));
  size_t mask = (size_t)key_len >= bsizeof(size_t) ? (size_t)-1 : ((size_t)1 << key_len) - 1;
  for(uint64_t i = start_pos; i < start_pos + bsizeof(uint64_t); ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    // Create mer m so that it will hash to position i
    m.randomize();
    m2 = m;
    m2.set_bits(0, ary.matrix().r(), (uint64_t)i);
    m.set_bits(0, ary.matrix().r(), ary.inverse_matrix().times(m2));

    // Add this one element to the hash
    ary.clear();
    bool         is_new      = false;
    size_t       id          = (size_t)-1;
    unsigned int carry_shift = 0;
    EXPECT_TRUE(ary.add(m, i, &carry_shift, &is_new, &id));
    EXPECT_TRUE(is_new);
    // Only expected to agree on the length of the key. Applies only
    // if key_len < lsize. The bits above key_len are pseudo-random
    EXPECT_EQ((size_t)i & mask, id & mask);

    // Every position but i in the hash should be empty
    uint64_t val;
    for(ssize_t j = -bsizeof(uint64_t); j <= (ssize_t)bsizeof(uint64_t); ++j) {
      SCOPED_TRACE(::testing::Message() << "j:" << j);
      val = (uint64_t)-1;
      size_t jd = (start_pos + j) & ary_size_mask;
      ASSERT_EQ(jd == id, ary.get_key_val_at_id(jd, get_mer, val) == large_array::FILLED);
      if(jd == id) {
        ASSERT_EQ(m2, get_mer);
        ASSERT_EQ((uint64_t)jd, val);
      }
    }
  }
}

TEST_P(HashArray, Collisions) {
  static const int nb_collisions = 4;
  std::vector<mer_dna> mers(nb_collisions);
  std::vector<mer_dna> mers2(nb_collisions);
  std::map<mer_dna, uint64_t> map;
  ASSERT_EQ((unsigned int)key_len / 2, mer_dna::k());

  SCOPED_TRACE(::testing::Message() << "key_len:" << key_len << " val_len:" << val_len << " reprobe:" << reprobe_limit);

  mers[0].polyA(); mers2[0].polyA();
  mers[1].polyC(); mers2[1].polyC();
  mers[2].polyG(); mers2[2].polyG();
  mers[3].polyT(); mers2[3].polyT();

  size_t start_pos = random() % (ary_size - bsizeof(uint64_t));
  for(uint64_t i = start_pos; i < start_pos + bsizeof(uint64_t); ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    ary.clear();
    map.clear();

    // Add mers that it will all hash to position i
    for(int j = 0; j < nb_collisions; ++j) {
      mers2[j].set_bits(0, ary.matrix().r(), (uint64_t)i);
      mers[j].set_bits(0, ary.matrix().r(), ary.inverse_matrix().times(mers2[j]));
      ary.add(mers[j], j);
      map[mers[j]] += j;
    }

    lazy_iterator it    = ary.iterator_all<lazy_iterator>();
    size_t        count = 0;
    while(it.next()) {
      SCOPED_TRACE(::testing::Message() << "it.key():" << it.key());
      ASSERT_FALSE(map.end() == map.find(it.key()));
      EXPECT_EQ(map[it.key()], it.val());
      ++count;
    }
    EXPECT_EQ(map.size(), count);
  }
}

struct arrays_type {
  large_array array;
  mer_map     map;
  arrays_type(size_t size, uint16_t key_len, uint16_t val_len, uint16_t reprobe_limit) :
    array(size, key_len, val_len, reprobe_limit), map()
  { }
};

typedef std::unique_ptr<arrays_type> arrays_ptr;
arrays_ptr fill_array(size_t nb_elts, size_t size, int key_len, int val_len, int reprobe_limit) {
  arrays_ptr arrays(new arrays_type(size, key_len, val_len, reprobe_limit));
  large_array& ary = arrays->array;
  mer_map&     map = arrays->map;

  mer_dna mer;
  for(int i = 0; i < nb_elts; ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    mer.randomize();
    map[mer] += i;
    // If get false, hash array filled up: double size
    bool res = ary.add(mer, i);
    if(!res) {
      // std::cerr << "Double size (" << size << " -> " << (2 * size) << ") nb_elts:" << nb_elts
      //           << " key_len:" << key_len << " val_len:" << val_len
      //           << " mer:" << mer << std::endl;
      // return std::make_pair(std::move(ary), std::move(map));
      return fill_array(nb_elts, 2 * size, key_len, val_len, reprobe_limit);
    }
  }
  return arrays;
}

TEST_P(HashArray, Iterator) {
  static const int nb_elts = 1 << (ary_lsize - 1 - (val_len == 1));
  SCOPED_TRACE(::testing::Message() << "key_len:" << key_len << " val_len:" << val_len << " reprobe:" << reprobe_limit);

  arrays_ptr   res = fill_array(nb_elts, ary_size, key_len, val_len, reprobe_limit);
  large_array& ary = res->array;
  mer_map &    map = res->map;

  eager_iterator it     = ary.iterator_all<eager_iterator>();
  lazy_iterator  lit    = ary.iterator_all<lazy_iterator>();
  stl_iterator   stl_it = ary.iterator_all<stl_iterator>();
  int count = 0;
  for( ; it.next(); ++stl_it) {
    ASSERT_TRUE(lit.next());
    ASSERT_NE(ary.end(), stl_it);
    mer_map::const_iterator mit = map.find(it.key());
    SCOPED_TRACE(::testing::Message() << "key:" << it.key());
    ASSERT_NE(map.end(), mit);
    EXPECT_EQ(mit->first, it.key());
    EXPECT_EQ(mit->second, it.val());
    EXPECT_EQ(mit->first, lit.key());
    EXPECT_EQ(mit->second, lit.val());
    EXPECT_EQ(mit->first, stl_it->first);
    EXPECT_EQ(mit->second, stl_it->second);
    EXPECT_EQ(it.id(), lit.id());
    EXPECT_EQ(it.id(), stl_it.id());
    ++count;
  }
  EXPECT_FALSE(lit.next());
  EXPECT_EQ(ary.end(), stl_it);
  EXPECT_EQ(map.size(), (size_t)count);

  count               = 0;
  const int nb_slices = 1;
  for(int i = 0; i < nb_slices; ++i) {
    SCOPED_TRACE(::testing::Message() << "slice:" << i << " nb_slices:" << nb_slices);
    region_iterator rit = ary.iterator_slice<region_iterator>(i, nb_slices);
    while(rit.next()) {
      ASSERT_GE(rit.oid(), rit.start());
      ASSERT_LT(rit.oid(), rit.end());
      mer_map::const_iterator mit = map.find(rit.key());
      ASSERT_NE(map.end(), mit);
      EXPECT_EQ(mit->first, rit.key());
      EXPECT_EQ(mit->second, rit.val());
      ++count;
    }
  }
  EXPECT_EQ(map.size(), (size_t)count);

  int i = 0;
  for(mer_map::const_iterator it = map.begin(); it != map.end(); ++it, ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i << " key:" << it->first);
    uint64_t val;
    size_t   id;
    EXPECT_TRUE(ary.get_key_id(it->first, &id));
    ASSERT_TRUE(ary.get_val_for_key(it->first, &val));
    EXPECT_EQ(it->second, val);
  }
}

TEST_P(HashArray, LargeValue) {
  mer_dna mer;
  mer.randomize();
  ary.add(mer, numeric_limits<uint64_t>::max());

  uint64_t val = 0;
  ASSERT_TRUE(ary.get_val_for_key(mer, &val));
  ASSERT_EQ(numeric_limits<uint64_t>::max(), val);
}

INSTANTIATE_TEST_CASE_P(HashArrayTest, HashArray, ::testing::Combine(::testing::Range(8, 4 * 64, 2), // Key lengths
                                                                     ::testing::Range(1, 10),    // Val lengths
                                                                     ::testing::Range(6, 8)      // Reprobe lengths
                                                                     ));

TEST(Hash, Set) {
  static const int lsize = 16;
  static const int size = 1 << lsize;
  static const int nb_elts = 2 * size / 3;

  large_array ary(size, 100, 0, 126);
  mer_set     set;
  mer_dna::k(50);
  mer_dna     mer;

  for(int i = 0; i < nb_elts; ++i) {
    mer.randomize();
    bool   is_new;
    size_t id;
    ASSERT_TRUE(ary.set(mer, &is_new, &id));
    ASSERT_EQ(set.insert(mer).second, is_new);
  }

  mer_dna tmp_mer;
  for(mer_set::const_iterator it = set.begin(); it != set.end(); ++it) {
    SCOPED_TRACE(::testing::Message() << "key:" << *it);
    size_t   id;
    EXPECT_TRUE(ary.get_key_id(*it, &id, tmp_mer));
  }

  for(int i = 0; i < nb_elts; ++i) {
    mer.randomize();
    size_t id;
    EXPECT_EQ(set.find(mer) != set.end(), ary.get_key_id(mer, &id));
  }
}

TEST(Hash, Update) {
  static const int lsize = 16;
  static const int size = 1 << lsize;
  static const int nb_elts = 2 * size / 3;

  large_array ary(size, 100, 4, 126);
  mer_map     in_ary;
  mer_dna::k(50);
  mer_dna     mer;

  for(int i = 0; i < nb_elts; ++i) {
    mer.randomize();
    bool is_new;
    size_t id;
    ASSERT_TRUE(ary.set(mer, &is_new, &id));
    auto res = in_ary.insert(std::make_pair(mer, (uint64_t)0));
    ASSERT_EQ(res.second, is_new);
  }

  for(auto it = in_ary.begin(); it != in_ary.end(); ++it) {
    uint64_t val = random_bits(4);
    EXPECT_TRUE(ary.update_add(it->first, val));
    it->second = val;
  }

  for(int i = 0; i < nb_elts; ++i) {
    mer.randomize();
    uint64_t val = random_bits(4);
    auto it = in_ary.find(mer);
    if(it == in_ary.end()) {
      EXPECT_FALSE(ary.update_add(mer, val));
    } else {
      it->second += val;
      EXPECT_TRUE(ary.update_add(mer, val));
    }
  }

  lazy_iterator it = ary.iterator_all<lazy_iterator>();
  size_t count = 0;
  while(it.next()) {
    ASSERT_NE(in_ary.end(), in_ary.find(it.key()));
    EXPECT_EQ(in_ary[it.key()], it.val());
    ++count;
  }
  EXPECT_EQ(in_ary.size(), count);
}

TEST(Hash, Info) {
  for(int iteration = 0; iteration < 100; ++iteration) {
    size_t                  mem     = random_bits(48);
    uint16_t                key_len = random_bits(7) + 1;
    uint16_t                val_len = random_bits(4) + 1;
    large_array::usage_info info(key_len, val_len, 126);

    SCOPED_TRACE(::testing::Message() << "iteration:" << iteration << " mem:" << mem
                 << " key_len:" << key_len << " val_len:" << val_len);
    uint16_t size_bits = info.size_bits(mem);
    uint16_t size2_bits = info.size_bits_linear(mem);
    ASSERT_EQ(size2_bits, size_bits);
    ASSERT_LE(info.mem((size_t)1 << size_bits), mem);
    ASSERT_GT(info.mem((size_t)1 << (size_bits + 1)), mem);
  }
}
}
