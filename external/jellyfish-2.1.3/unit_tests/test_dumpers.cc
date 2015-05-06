#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/file_header.hpp>
#include <jellyfish/binary_dumper.hpp>
#include <jellyfish/text_dumper.hpp>
#include <jellyfish/mapped_file.hpp>

namespace {
using jellyfish::mer_dna;
using jellyfish::file_header;
typedef jellyfish::cooperative::hash_counter<mer_dna> hash_counter;
typedef hash_counter::array::eager_iterator iterator;

struct binary {
  typedef jellyfish::binary_dumper<hash_counter::array> dumper;
  typedef jellyfish::binary_reader<mer_dna, uint64_t> reader;
  typedef jellyfish::binary_query_base<mer_dna, uint64_t> query;
};
struct text {
  typedef jellyfish::text_dumper<hash_counter::array> dumper;
  typedef jellyfish::text_reader<mer_dna, uint64_t> reader;
};


TEST(Dumper, IO) {
  static const int   mer_len          = 50;
  static const int   hash_size        = 5000;
  static const int   nb               = hash_size;
  static const int   hash_val_len     = 5;
  static const int   dump_counter_len = 1;
  static const char* file_binary      = "./binary_dumper";
  static const char* file_text        = "./text_dumper";

  file_unlink bf(file_binary);
  file_unlink tf(file_text);

  mer_dna::k(mer_len);
  hash_counter hash(hash_size, mer_len * 2, 5 /* val len */, 1 /* nb threads */);

  mer_dna m;
  for(int i = 0; i < nb; i++) {
    m.randomize();
    const uint64_t rval = random_bits(9);
    hash.add(m, rval);
    uint64_t val;
    ASSERT_TRUE(hash.ary()->get_val_for_key(m, &val));
    EXPECT_EQ(rval, val);
  }

  // Dump without zeroing to check dumped content against in memory hash
  {
    file_header bh;
    bh.fill_standard();
    bh.update_from_ary(*hash.ary());
    binary::dumper bd(dump_counter_len, mer_len * 2, 4, file_binary, &bh);
    bd.one_file(true);
    bd.zero_array(false);
    bd.dump(hash.ary());

    file_header th;
    th.fill_standard();
    th.update_from_ary(*hash.ary());
    text::dumper td(4, file_text, &th);
    td.one_file(true);
    td.zero_array(false);
    td.dump(hash.ary());
  }

  // Check dumped content
  {
    file_header bh;
    std::ifstream bis(file_binary);
    bh.read(bis);
    EXPECT_STREQ(binary::dumper::format, bh.format().c_str());
    EXPECT_EQ(dump_counter_len, bh.counter_len());
    binary::reader br(bis, &bh);

    jellyfish::mapped_file binary_map(file_binary);
    binary::query bq(binary_map.base() + bh.offset(), bh.key_len(), bh.counter_len(), bh.matrix(),
                     bh.size() - 1, binary_map.length() - bh.offset());

    file_header th;
    std::ifstream tis(file_text);
    th.read(tis);
    EXPECT_STREQ(binary::dumper::format, bh.format().c_str());
    text::reader tr(tis, &th);

    const uint64_t max_val = ((uint64_t)1 << (8 * dump_counter_len)) - 1;
    int bcount = 0, tcount = 0, qcount = 0;
    mer_dna tmp_key;
    while(br.next()) {
      uint64_t val = 0;
      size_t   id  = 0;
      bool present = hash.ary()->get_val_for_key(br.key(), &val, tmp_key, &id);
      EXPECT_TRUE(present);
      if(present) {
        EXPECT_EQ(std::min(max_val, val), br.val());
        ++bcount;
      }

      EXPECT_TRUE(tr.next());
      present = hash.ary()->get_val_for_key(tr.key(), &val);
      EXPECT_TRUE(present);
      if(present) {
        EXPECT_EQ(val, tr.val());
        ++tcount;
      }

      uint64_t query_val;
      uint64_t query_id;
      present = bq.val_id(br.key(), &query_val, &query_id);
      EXPECT_TRUE(present);
      if(present) {
        EXPECT_EQ(std::min(max_val, val), query_val);
        // EXPECT_EQ(id, query_id);
        ++qcount;
      }
    }
    EXPECT_EQ(nb, bcount);
    EXPECT_EQ(nb, tcount);
    EXPECT_EQ(nb, qcount);
  }

  // Dump with zeroing and check hash is empty
  {
    file_header bh;
    bh.fill_standard();
    bh.update_from_ary(*hash.ary());
    binary::dumper bd(dump_counter_len, mer_len * 2, 4, file_binary, &bh);
    bd.one_file(true);
    bd.zero_array(true);
    bd.dump(hash.ary());
  }
  {
    iterator it = hash.ary()->eager_slice(0, 1);
    uint64_t count = 0;
    while(it.next()) ++count;
    EXPECT_EQ((uint64_t)0, count);
  }
}

} // namespace {
