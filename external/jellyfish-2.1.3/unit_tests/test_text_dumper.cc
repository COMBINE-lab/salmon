#include <map>
#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <jellyfish/text_dumper.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/thread_exec.hpp>
#include <jellyfish/mer_dna.hpp>

namespace {
using jellyfish::mer_dna;
typedef jellyfish::cooperative::hash_counter<mer_dna> hash_counter;
typedef hash_counter::array::lazy_iterator lazy_iterator;
typedef hash_counter::array::eager_iterator eager_iterator;
typedef hash_counter::array::region_iterator region_iterator;
typedef jellyfish::text_dumper<hash_counter::array> dumper;

class hash_adder : public jellyfish::thread_exec {
  typedef std::map<mer_dna, uint64_t> map;
  typedef std::vector<map>            maps;

  hash_counter& hash_;
  int           nb_;
  maps          check_;

public:
  hash_adder(hash_counter& hash, int nb, int nb_threads) :
    hash_(hash),
    nb_(nb),
    check_(nb_threads)
  { }
  void start(int id) {
    mer_dna m;
    map&    my_map = check_[id];

    for(int i = 0; i < nb_; ++i) {
      m.randomize();
      hash_.add(m, 1);
      ++my_map[m];
    }

    hash_.done();
  }

  void maps_consolidate() {
    maps::iterator first = check_.begin();
    maps::const_iterator it = first;
    for(++it; it < check_.end(); ++it)
      for(map::const_iterator vit = it->begin(); vit != it->end(); ++vit)
        (*first)[vit->first] += vit->second;
  }

  void hash_check() {
    map& m = check_.front();
    size_t id;
    for(map::const_iterator vit = m.begin(); vit != m.end(); ++vit)
      EXPECT_TRUE(hash_.ary()->get_key_id(vit->first, &id));
  }

  uint64_t operator[](const mer_dna& m) const {
    map::const_iterator it = check_.front().find(m);
    if(it == check_.front().end())
      return 0;
    else
      return it->second;
  }

  size_t size() const { return check_.front().size(); }
 };

template<typename T>
struct inc_on_destructor {
  T* const x;
  inc_on_destructor(T* nb) : x(nb) { }
  ~inc_on_destructor() { ++*x; }
  operator T() const { return *x; }
};

TEST(TextDumper, Random) {
  static int iteration_count_ = 0;
  ++iteration_count_;

  static const int    mer_len    = 35;
  static const int    nb_threads = 5;
  static const int    nb         = 5000;
  static const size_t init_size  = 5000;
  static const char*  file       = "./text_dumper";
  file_unlink f(file);

  mer_dna::k(mer_len);
  hash_counter hash(init_size, mer_len * 2, 5, nb_threads, 15);
  hash_adder adder(hash, nb, nb_threads);
  adder.exec_join(nb_threads);
  adder.maps_consolidate();
  adder.hash_check();

  {
    dumper d(nb_threads, file);
    d.one_file(true);
    d.dump(hash.ary());
  }

  std::ifstream fd(file);
  std::string line;
  mer_dna m;
  uint64_t count = 0;
  while(std::getline(fd, line)) {
    m = line.substr(0, mer_dna::k());
    EXPECT_EQ(adder[m], (uint64_t)atol(line.substr(mer_dna::k() + 1).c_str()));
    ++count;
  }
  EXPECT_EQ(adder.size(), count);

  //  unlink(file);
}
} // namespace {
