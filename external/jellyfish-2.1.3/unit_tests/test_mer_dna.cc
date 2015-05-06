/* SuperRead pipeline
 * Copyright (C) 2012  Genome group at University of Maryland.
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#include <stdio.h>
#include <map>
#include <iostream>
#include <sstream>
#include <gtest/gtest.h>
#include <jellyfish/mer_dna.hpp>
#include <unit_tests/test_main.hpp>

namespace {
using namespace jellyfish;

const std::string short_mer("ACGTTGGCCAA");
const std::string mid_mer("AAGACTTAGAATCAGCTAAGGTAACTTACGTAGAATATAGA");
const std::string long_mer("AAGACTTAGAATCAGCTAAGGTAACTTACGTAGAATATAG"
                           "AAGACGGCGCCTATCCCAACCTCCATTTGCTGACGCCCTC"
                           "GTAACCGTGCTGCGAGGTTACTCTATACTGA");
const std::string huge_mer("CAGGACACAATACGTCGATCCAATCCCCGACGTGAGAGTT"
                           "TAACGCTAATCTTGATCCATTACAAGTATAGATATTTCGG"
                           "GCGCCACGGGGAAACTTGCCCTATGTTCGAGTCCGCCACC"
                           "GCGGCACCAGCCTTTGTGTGACGGCCACAAAGGGTAAAAG"
                           "ATGTTGTCTCGCGCCGGCGTGCGGTCTTTACCAGATCCTT"
                           "GAGCGGCCTAGAAAGTTCGTACCAGTTCGACTGAAACAAG"
                           "ACAACGGATCGCCCACGGTCATCACACGCCCACGCACGGG"
                           "GTCGGGTTGGTATATTCAACCTCGAGTTAAACGT");
const std::string* test_mers[4] = {
  &short_mer, &mid_mer, &long_mer, &huge_mer
};

TEST(MerDNASimple, InitSize64) {
  struct sinfo {
    unsigned int k;
    unsigned int nb_words;
    unsigned int nb_msb;
    uint64_t     msw;
    unsigned int lshift;
  };
  sinfo ary[5] = {
    { 5, 1, 10, (uint64_t)0x3ff, 8 }, { 32, 1, 64, (uint64_t)-1, 62 },
    { 33, 2, 2, (uint64_t)0x3, 0 }, { 64, 2, 64, (uint64_t)-1, 62 },
    { 65, 3, 2, (uint64_t)0x3, 0 }
  };
  typedef mer_dna_ns::mer_base_dynamic<uint64_t> mer64;
  for(size_t i = 0; i < sizeof(ary) / sizeof(sinfo); ++i) {
    mer64 m(ary[i].k);
    EXPECT_EQ(ary[i].k, m.k());
    EXPECT_EQ(ary[i].nb_words, m.nb_words());
    EXPECT_EQ(ary[i].nb_msb, m.nb_msb());
    EXPECT_EQ(ary[i].lshift, m.lshift());
  }
}

#ifdef HAVE_INT128
TEST(MerDNASimple, InitSize128) {
  struct sinfo {
    unsigned int      k;
    unsigned int      nb_words;
    unsigned int      nb_msb;
    unsigned __int128 msw;
    unsigned int      lshift;
  };
  sinfo ary[5] = {
    { 5, 1, 10, (unsigned __int128)0x3ff, 8 },
    { 32, 1, 64, (unsigned __int128)0xffffffffffffffffUL, 62 },
    { 33, 1, 66, (unsigned __int128)0xffffffffffffffffUL | ((unsigned __int128)0x3 << 64), 64 },
    { 64, 1, 128, (unsigned __int128)-1, 126 },
    { 65, 2, 2, (unsigned __int128)0x3, 0 }
  };
  typedef mer_dna_ns::mer_base_dynamic<unsigned __int128> mer128;
  for(size_t i = 0; i < sizeof(ary) / sizeof(sinfo); ++i) {
    mer128 m(ary[i].k);
    EXPECT_EQ(ary[i].k, m.k());
    EXPECT_EQ(ary[i].nb_words, m.nb_words());
    EXPECT_EQ(ary[i].nb_msb, m.nb_msb());
    EXPECT_EQ(ary[i].lshift, m.lshift());
  }
}
#endif

TEST(MerDNASimple, Codes) {
  EXPECT_EQ(mer_dna::CODE_A, mer_dna::code('A'));
  EXPECT_EQ(mer_dna::CODE_A, mer_dna::code('a'));
  EXPECT_EQ(mer_dna::CODE_C, mer_dna::code('C'));
  EXPECT_EQ(mer_dna::CODE_C, mer_dna::code('c'));
  EXPECT_EQ(mer_dna::CODE_G, mer_dna::code('G'));
  EXPECT_EQ(mer_dna::CODE_G, mer_dna::code('g'));
  EXPECT_EQ(mer_dna::CODE_T, mer_dna::code('T'));
  EXPECT_EQ(mer_dna::CODE_T, mer_dna::code('t'));
  EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_A));
  EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_C));
  EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_G));
  EXPECT_FALSE(mer_dna::not_dna(mer_dna::CODE_T));

  for(int c = 0; c < 256; ++c) {
    switch((char)c) {
    case 'A': case 'a':
    case 'C': case 'c':
    case 'G': case 'g':
    case 'T': case 't':
      EXPECT_FALSE(mer_dna::not_dna(mer_dna::code(c)));
      break;
    default:
      EXPECT_TRUE(mer_dna::not_dna(mer_dna::code(c)));
      break;
    }
  }
}

TEST(MerDNASimple, SetBits) {
  mer_dna::k(100);
  static const int pattern_len = 10;
  uint64_t pattern = ::random_bits(pattern_len); // Create a random pattern
  mer_dna even_pattern, odd_pattern;             // And the corresponding mers
  for(int i = pattern_len - 2; i >= 0; i -= 2) {
    even_pattern.shift_left((int)((pattern >> i) & 0x3));
    odd_pattern.shift_left((int)((pattern >> (i + 1)) & 0x3));
  }
  odd_pattern.shift_left((int)((pattern << 1) & 0x3));

  mer_dna mer;
  for(int i = 0; i <= (int)(mer_dna::k() - pattern_len / 2); ++i, even_pattern.shift_left('A'), odd_pattern.shift_left('A')) {
    // Even
    mer.polyA();
    EXPECT_EQ(std::string(mer_dna::k(), 'A'), mer.to_str());
    mer.set_bits(2 * i, pattern_len, pattern);
    EXPECT_EQ(pattern, mer.get_bits(2 * i, pattern_len));
    EXPECT_EQ(even_pattern, mer);
    EXPECT_EQ(even_pattern.to_str(), mer.to_str());
    // Odd
    mer.polyA();
    mer.set_bits(2 * i + 1, pattern_len, pattern);
    if(i < (int)(mer_dna::k() - pattern_len / 2))
      EXPECT_EQ(pattern, mer.get_bits(2 * i + 1, pattern_len));
    else // On the largest value of i, one bit may have fallen off the end of the mer
      EXPECT_EQ(pattern & (((uint64_t)1 << (pattern_len - 1)) - 1), mer.get_bits(2 * i + 1, pattern_len));
    EXPECT_EQ(odd_pattern, mer);
    EXPECT_EQ(odd_pattern.to_str(), mer.to_str());
  }
}

TEST(MerDNASimple, Shifts) {
  for(int i = 1; i < 100; ++i) {
    mer_dna::k(i);
    mer_dna m;

    m.randomize();
    const int c = ::random_bits(2);
    mer_dna rm(m);
    rm.shift_right(c);
    mer_dna lm(m);
    lm.shift_left(c);

    EXPECT_EQ(c, rm.base(mer_dna::k() - 1).code());
    EXPECT_EQ(c, lm.base(0).code());
    for(unsigned int j = 0; j < mer_dna::k() - 1; ++j) {
      EXPECT_EQ(m.base(j + 1).code(), rm.base(j).code());
      EXPECT_EQ(m.base(j).code(), lm.base(j + 1).code());
    }
  }
} // MerDNASimple.Shifts


bool simple_homolymer_test(const mer_dna& m) {
  mer_dna cm(m);
  cm.shift_right(m.base(0).code());
  return cm == m;
}

TEST(MerDNASimple, Homopolymer) {
  for(int i = 1; i < 256; ++i) {
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    mer_dna::k(i);
    mer_dna m;

    for(int j = 0; j < 10; ++j) {
      m.randomize();
      EXPECT_EQ(simple_homolymer_test(m), m.is_homopolymer());
    }

    m.polyA();
    EXPECT_TRUE(simple_homolymer_test(m));
    EXPECT_TRUE(m.is_homopolymer());
    if(i > 1) { // i == 1 all mers are homopolymers by definition
      m.base(::random_bits(5) % i) = 'T';
      if(simple_homolymer_test(m) || m.is_homopolymer())
        std::cerr << m << "\n";
      EXPECT_FALSE(simple_homolymer_test(m));
      EXPECT_FALSE(m.is_homopolymer());
    }

    m.polyC();
    EXPECT_TRUE(simple_homolymer_test(m));
    EXPECT_TRUE(m.is_homopolymer());
    m.polyG();
    EXPECT_TRUE(simple_homolymer_test(m));
    EXPECT_TRUE(m.is_homopolymer());
    m.polyT();
    EXPECT_TRUE(simple_homolymer_test(m));
    EXPECT_TRUE(m.is_homopolymer());
  }
} // MerDNASimple.Homopolymer

TEST(MerDNASimple, Comparators) {
  mer_dna::k(151);
  mer_dna ma, mc, mg, mt;
  ma.polyA();
  mc.polyC();
  mg.polyG();
  mt.polyT();

  mer_dna cma(ma), cmc(mc), cmg(mg), cmt(mt);

  ASSERT_TRUE(ma < mc); ASSERT_FALSE(mc < ma);
  ASSERT_TRUE(ma < mg); ASSERT_FALSE(mg < ma);
  ASSERT_TRUE(ma < mt); ASSERT_FALSE(mt < ma);
  ASSERT_TRUE(mc < mg); ASSERT_FALSE(mg < mc);
  ASSERT_TRUE(mc < mt); ASSERT_FALSE(mt < mc);
  ASSERT_TRUE(mg < mt); ASSERT_FALSE(mt < mg);

  ASSERT_FALSE(ma < ma);
  ASSERT_FALSE(mc < mc);
  ASSERT_FALSE(mg < mg);
  ASSERT_FALSE(mt < mt);

  std::map<mer_dna, int> map;
  EXPECT_EQ(1, ++map[ma]);
  EXPECT_EQ(1, ++map[mc]);
  EXPECT_EQ(1, ++map[mg]);
  EXPECT_EQ(1, ++map[mt]);
  EXPECT_EQ(2, ++map[cma]);
  EXPECT_EQ(2, ++map[cmc]);
  EXPECT_EQ(2, ++map[cmg]);
  EXPECT_EQ(2, ++map[cmt]);

  mer_dna m1, m2;
  for(int i = 0; i < 1000; ++i) {
    m1.randomize();
    m2.randomize();
    SCOPED_TRACE(::testing::Message() << "m1:" << m1 << " m2:" << m2);
    ASSERT_FALSE(m1 == m2); // Very small probability to fail (k == 151)

    EXPECT_TRUE(m1 == m1);
    EXPECT_FALSE(m1 < m1);
    EXPECT_FALSE(m1 > m1);
    EXPECT_TRUE(m1 <= m1);
    EXPECT_TRUE(m1 >= m1);

    EXPECT_TRUE(m2 == m2);
    EXPECT_FALSE(m2 < m2);
    EXPECT_FALSE(m2 > m2);
    EXPECT_TRUE(m2 <= m2);
    EXPECT_TRUE(m2 >= m2);

    EXPECT_EQ(m1.to_str().compare(m2.to_str()) < 0, m1 < m2); // Comparison is lexicographic

    EXPECT_TRUE(m1 < m2 || m2 < m1);
    EXPECT_TRUE(m1 <= m2 || m2 <= m1);
    EXPECT_FALSE(m1 <= m2 && m2 <= m1);
    EXPECT_TRUE(m1 > m2 || m2 > m1);
    EXPECT_TRUE(m1 >= m2 || m2 >= m1);
    EXPECT_FALSE(m1 >= m2 && m2 >= m1);

    EXPECT_NE(m1 < m2, m1 >= m2);
    EXPECT_NE(m1 < m2, m1 > m2);
  }
}

TEST(MerDNASimple, IO) {
  std::stringstream buffer;

  for(int i = 0; i < 10000; ++i) {
    buffer.clear();
    SCOPED_TRACE(::testing::Message() << "i:" << i);
    mer_dna::k(::random_bits(9) + 1);
    mer_dna m1, m2;
    m1.randomize();
    buffer << m1;
    buffer >> m2;
    EXPECT_EQ(m1, m2);
  }
}

TEST(MerDNASimple, MultipleSize) {
  typedef jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1> mer_dna1;
  typedef jellyfish::mer_dna_ns::mer_base_static<uint64_t, 2> mer_dna2;

  mer_dna::k(10);
  mer_dna1::k(50);
  mer_dna2::k(100);
  EXPECT_EQ(10, mer_dna::k());
  EXPECT_EQ(0, mer_dna::class_index);
  EXPECT_EQ(50, mer_dna1::k());
  EXPECT_EQ(1, mer_dna1::class_index);
  EXPECT_EQ(100, mer_dna2::k());
  EXPECT_EQ(2, mer_dna2::class_index);
}

// Value Type Container class
template <typename T, int N>
class VTC {
public:
  typedef T Type;
  static const int test_id = N;
};
template <typename T, int N>
const int VTC<T, N>::test_id;

template<typename VT>
class MerDNA : public ::testing::Test {
public:
  typedef typename VT::Type Type;
  void SetUp() {
    Type::k(GetParam().size());
  }
  const std::string& GetParam() const {
    return *test_mers[VT::test_id];
  }
};
typedef ::testing::Types<VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 0>,
                         VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 1>,
                         VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 2>,
                         VTC<mer_dna_ns::mer_base_dynamic<uint64_t>, 3>,
#ifdef HAVE_INT128
                         VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 0>,
                         VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 1>,
                         VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 2>,
                         VTC<mer_dna_ns::mer_base_dynamic<unsigned __int128>, 3>,
#endif
                         VTC<mer_dna_ns::mer_base_static<uint32_t>, 3>,
                         VTC<mer_dna_ns::mer_base_static<uint64_t>, 0>,
                         VTC<mer_dna_ns::mer_base_static<uint64_t>, 1>,
                         VTC<mer_dna_ns::mer_base_static<uint64_t>, 2>,
                         VTC<mer_dna_ns::mer_base_static<uint64_t>, 3>,
#ifdef HAVE_INT128
                         VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 0>,
                         VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 1>,
                         VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 2>,
                         VTC<mer_dna_ns::mer_base_static<unsigned __int128>, 3>,
#endif
                         VTC<mer_dna_ns::mer_base_static<uint32_t>, 3>
                         > MerDNATypes;
TYPED_TEST_CASE(MerDNA, MerDNATypes);

TYPED_TEST(MerDNA, InitFromStr) {
  typename TypeParam::Type m(this->GetParam());
  EXPECT_EQ(this->GetParam().size(), m.k());
  EXPECT_EQ(this->GetParam(), m.to_str());
}

TYPED_TEST(MerDNA, ShiftLeft) {
  typename TypeParam::Type m(this->GetParam().size());
  m.polyA();
  int inserted = 0;
  for(std::string::const_iterator it = this->GetParam().begin(); it != this->GetParam().end(); ++it, ++inserted) {
    m.shift_left(*it);

    int check = inserted;
    for(std::string::const_iterator cit = this->GetParam().begin(); check >= 0; ++cit, --check)
      EXPECT_EQ(*cit, (char)m.base(check));
  }
  EXPECT_EQ(this->GetParam(), m.to_str());
}

TYPED_TEST(MerDNA, ShiftRight) {
  typename TypeParam::Type m(this->GetParam().size());
  m.polyA();
  int inserted = 0;
  for(std::string::const_reverse_iterator it = this->GetParam().rbegin(); it != this->GetParam().rend(); ++it, ++inserted) {
    m.shift_right(*it);

    int check = inserted;
    for(std::string::const_reverse_iterator cit = this->GetParam().rbegin(); check >= 0; ++cit, --check)
      EXPECT_EQ(*cit, (char)m.base(m.k() - 1 - check));
  }
  EXPECT_EQ(this->GetParam(), m.to_str());
}

TYPED_TEST(MerDNA, Equality) {
  typename TypeParam::Type m1(this->GetParam());
  typename TypeParam::Type m2(this->GetParam().size());

  char str[this->GetParam().size() + 1];
  str[this->GetParam().size()] = '\0';
  memset(str, 'A', this->GetParam().size());
  m2.polyA();
  EXPECT_STREQ(str, m2.to_str().c_str());
  memset(str, 'C', this->GetParam().size());
  m2.polyC();
  EXPECT_STREQ(str, m2.to_str().c_str());
  memset(str, 'G', this->GetParam().size());
  m2.polyG();
  EXPECT_STREQ(str, m2.to_str().c_str());
  memset(str, 'T', this->GetParam().size());
  m2.polyT();
  EXPECT_STREQ(str, m2.to_str().c_str());


  int i = 1;
  for(std::string::const_iterator it = this->GetParam().begin(); it < this->GetParam().end(); ++it, ++i) {
    sprintf(str + this->GetParam().size() - i, "%.*s", i, this->GetParam().c_str());
    typename TypeParam::Type m(str);
    EXPECT_STREQ(str, m.to_str().c_str());
    m2.shift_left(*it);
    EXPECT_EQ(m, m2);
  }
  EXPECT_TRUE(m1 == m2);
  EXPECT_FALSE(m1 != m2);
  EXPECT_EQ(m1.to_str(), m2.to_str());

  // typename TypeParam::Type m3(this->GetParam());
  // m3[0] = 0;
  // EXPECT_FALSE(m1 == m3);

  // typename TypeParam::Type m4(this->GetParam().size() + 1);
  // EXPECT_FALSE(m1 == m4);
  // typename TypeParam::Type m5(this->GetParam().size() - 1);
  // EXPECT_FALSE(m1 == m5);

}

TYPED_TEST(MerDNA, Copy) {
  typename TypeParam::Type m1(this->GetParam());
  typename TypeParam::Type m2(m1);
  typename TypeParam::Type m3(this->GetParam().size());
  m3 = m1;

  EXPECT_TRUE(m1 == m2);
  EXPECT_TRUE(m2 == m3);
  EXPECT_TRUE(m3 == m1);
  m1.shift_left('A');
  EXPECT_TRUE(!(m1 == m2));
  EXPECT_TRUE(!(m1 == m3));
}

TYPED_TEST(MerDNA, OperatorShift) {
  typename TypeParam::Type m(this->GetParam());
  std::ostringstream os;
  os << m;
  EXPECT_EQ(this->GetParam(), os.str());
}

TYPED_TEST(MerDNA, GetBits) {
  typename TypeParam::Type m(this->GetParam());
  for(unsigned int i = 0; i < 20; ++i) {
    long int start   = random() % (this->GetParam().size() - 1);
    long int max_len =
      std::min(this->GetParam().size() - start, 8 * sizeof(typename TypeParam::Type::base_type));
    long int len     = (random() % (max_len - 1)) + 1;

    // Get bits by right-shifting
    typename TypeParam::Type cm(m);
    for(unsigned int j = 1; j < start; j += 2)
      cm.shift_right(0); // Shift by 2 bits
    typename TypeParam::Type::base_type y = cm.word(0);
    if(start & 0x1)
      y >>= 1;
    y &= ((typename TypeParam::Type::base_type)1 << len) - 1;

    EXPECT_EQ(y, m.get_bits(start, len));
  }
}
TYPED_TEST(MerDNA, GetBases) {
  typename TypeParam::Type m(this->GetParam());

  for(std::string::const_reverse_iterator it = this->GetParam().rbegin(); it != this->GetParam().rend(); ++it)
    EXPECT_EQ(*it, (char)m.base(it - this->GetParam().rbegin()));

  const char bases[4] = { 'A', 'C', 'G', 'T' };
  for(const char* it = bases; it != bases + 4; ++it) {
    typename TypeParam::Type n(m);
    for(size_t j = 0; j < this->GetParam().size(); ++j)
      n.base(j) = *it;
    typename TypeParam::Type m_expected(std::string(this->GetParam().size(), *it));
    EXPECT_EQ(m_expected, n);
  }
}

char rc_base(char c) {
  switch(c) {
  case 'A': case 'a': return 'T';
  case 'C': case 'c': return 'G';
  case 'G': case 'g': return 'C';
  case 'T': case 't': return 'A';
  }
  return 'A'; // Should never be reached
}

TYPED_TEST(MerDNA, ReverseComplement) {
  typename TypeParam::Type m(this->GetParam());
  std::string rc(this->GetParam().size(), 'A');
  for(size_t i = 0; i < rc.size(); ++i)
    rc[i] = rc_base(this->GetParam()[this->GetParam().size() - 1 - i]);
  EXPECT_EQ(this->GetParam().size(), m.k());
  m.reverse_complement();
  EXPECT_EQ(rc, m.to_str());
  typename TypeParam::Type rm(rc);
  EXPECT_EQ(rm, m);
  EXPECT_EQ(m, m.get_reverse_complement().get_reverse_complement());
}

TYPED_TEST(MerDNA, Canonical) {
  typename TypeParam::Type m(this->GetParam());
  typename TypeParam::Type canonical = m.get_canonical();

  EXPECT_FALSE(m < canonical);
  EXPECT_TRUE(canonical <= m);
  EXPECT_TRUE(canonical == m || canonical == m.get_reverse_complement());
  m.canonicalize();
  EXPECT_EQ(canonical, m.get_canonical());
}

} // namespace {
