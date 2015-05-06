#include <gtest/gtest.h>
#include <jellyfish/cooperative_pool2.hpp>
#include <jellyfish/thread_exec.hpp>

namespace {
// Generate all numbers in [0, producers * max)
class sequence : public jellyfish::cooperative_pool2<sequence, int> {
  typedef jellyfish::cooperative_pool2<sequence, int> super;

  const uint32_t   max_;
  std::vector<int> cur_;
  const uint32_t   producers_;

public:
  std::vector<int> check_;
  sequence(uint32_t producers, uint32_t threads, uint32_t max) :
    super(producers, 3 * threads),
    max_(max),
    cur_(producers, 0),
    producers_(producers),
    check_(max * producers, 0)
  { }

  bool produce(uint32_t i, int& e) {
    assert(i < producers_);
    int& cur = cur_[i];
    if(cur < max_) {
      e = i * max_ + cur++;
      __sync_add_and_fetch(&check_[e], 1);
      return false;
    }
    return true;
  }
};

class list_ints : public jellyfish::thread_exec {
  sequence         seq_;
  std::vector<int> check_;
public:
  list_ints(uint32_t producers, uint32_t threads, uint32_t max) : seq_(producers, threads, max), check_(producers * max, 0) { }
  virtual void start(int i) {
    while(true) {
      sequence::job j(seq_);
      if(j.is_empty())
        break;
      ++check_[*j];
    }
  }

  bool check() const {
    for(auto it = check_.cbegin(); it != check_.cend(); ++it)
      if(*it != 1)
        return false;
    return true;
  }

  bool check_print() const {
    if(check())
      return true;

    for(auto it = check_.cbegin(); it != check_.cend(); ++it)
      printf("%d", *it);
    printf("\n----------\n");
    for(auto it = seq_.check_.cbegin(); it != seq_.check_.cend(); ++it)
      printf("%d", *it);
    printf("\n");
    return false;
  }
};

class CooperativePoolTest : public ::testing::TestWithParam<uint32_t> {
public:
  static const uint32_t nb_threads = 10;
  CooperativePoolTest() : workers(GetParam(), nb_threads, 1000) { }

protected:
  list_ints workers;
};

TEST_P(CooperativePoolTest, Ints) {
  workers.exec_join(nb_threads);
  EXPECT_TRUE(workers.check_print());
}

INSTANTIATE_TEST_CASE_P(CooperativePool,
                        CooperativePoolTest,
                        ::testing::Range((uint32_t)1, CooperativePoolTest::nb_threads + 1));

} // namespace {
