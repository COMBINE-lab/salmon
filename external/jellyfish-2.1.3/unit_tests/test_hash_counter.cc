#include <gtest/gtest.h>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/mer_dna.hpp>
#include <map>
#include <vector>
#include <limits>

namespace {
using jellyfish::thread_exec;
using jellyfish::mer_dna;
typedef jellyfish::cooperative::hash_counter<mer_dna> hash_counter;
typedef hash_counter::array::lazy_iterator lazy_iterator;

enum OPERATION { ADD, SET };

class hash_adder : public thread_exec {
  typedef std::map<mer_dna, uint64_t> map;
  typedef std::vector<map>            maps;

  hash_counter& hash_;
  int           nb_;
  maps          check_;
  OPERATION     op_;

public:
  hash_adder(hash_counter& hash, int nb, int nb_threads, OPERATION op) :
    hash_(hash),
    nb_(nb),
    check_(nb_threads),
    op_(op)
  { }
  void start(int id) {
    mer_dna m;
    map&    my_map = check_[id];

    for(int i = 0; i < nb_; ++i) {
      m.randomize();
      switch(op_) {
      case ADD:
        hash_.add(m, std::numeric_limits<uint64_t>::max());
        break;
      case SET:
        hash_.set(m);
        break;
      }

      my_map[m] = std::numeric_limits<uint64_t>::max();
    }

    hash_.done();
  }

  uint64_t val(const mer_dna& m) {
    uint64_t res = 0;
    for(maps::const_iterator it = check_.begin(); it < check_.end(); ++it) {
      map::const_iterator vit = (*it).find(m);
      if(vit != it->end())
        res += vit->second;
    }
    return res;
  }
};

TEST(HashCounterCooperative, SizeDouble) {
  static const int    mer_len    = 35;
  static const int    nb_threads = 5;
  static const int    nb         = 200;
  static const size_t init_size  = 128;
  mer_dna::k(mer_len);

  {
    hash_counter hash(init_size, mer_len * 2, 5, nb_threads);
    EXPECT_TRUE(hash.do_size_doubling());
    EXPECT_EQ(mer_len * 2, hash.key_len());
    EXPECT_EQ(5, hash.val_len());

    hash_adder adder(hash, nb, nb_threads, ADD);
    adder.exec_join(nb_threads);

    lazy_iterator it = hash.ary()->iterator_all<lazy_iterator>();
    while(it.next())
      EXPECT_EQ(adder.val(it.key()), it.val());
    EXPECT_LT((size_t)(nb_threads * nb), hash.size());
  }

  {
    hash_counter hash(init_size, mer_len * 2, 0, nb_threads);
    EXPECT_TRUE(hash.do_size_doubling());
    EXPECT_EQ(mer_len * 2, hash.key_len());
    EXPECT_EQ(0, hash.val_len());

    hash_adder adder(hash, nb, nb_threads, SET);
    adder.exec_join(nb_threads);

    lazy_iterator it = hash.ary()->iterator_all<lazy_iterator>();
    while(it.next()) {
      SCOPED_TRACE(::testing::Message() << "mer:" << it.key());
      EXPECT_EQ(0, it.val());
    }
    EXPECT_LT((size_t)(nb_threads * nb), hash.size());
  }
}
} // namespace {
