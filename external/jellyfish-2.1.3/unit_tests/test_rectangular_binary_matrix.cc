#include <gtest/gtest.h>
#include <unit_tests/test_main.hpp>
#include <stdexcept>
#include <stdlib.h>
#include <jellyfish/rectangular_binary_matrix.hpp>
#include <jellyfish/mer_dna.hpp>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_INT128
#include <jellyfish/int128.hpp>
#endif

namespace {
using jellyfish::RectangularBinaryMatrix;
using jellyfish::mer_dna;

static bool allocate_matrix(unsigned int r, unsigned c) {
  RectangularBinaryMatrix m(r, c);
  return m.is_zero();
}

TEST(RectangularBinaryMatrix, InitSizes) {
  RectangularBinaryMatrix m(5, 60);

  EXPECT_EQ(5u, m.r());
  EXPECT_EQ(60u, m.c());
  EXPECT_TRUE(m.is_zero());

  EXPECT_THROW(allocate_matrix(100, 100), std::out_of_range);
  EXPECT_THROW(allocate_matrix(0, 100), std::out_of_range);
  EXPECT_THROW(allocate_matrix(10, 0), std::out_of_range);
}

TEST(RectangularBinaryMatrix, Copy) {
  RectangularBinaryMatrix m1(5, 55);
  m1.randomize(random_bits);

  RectangularBinaryMatrix m2(m1);
  RectangularBinaryMatrix m3(6, 66);
  RectangularBinaryMatrix m4(5, 55);

  EXPECT_TRUE(!m1.is_zero());
  EXPECT_TRUE(m1 == m2);
  EXPECT_TRUE(!(m1 == m3));
  EXPECT_TRUE(!(m1 == m4));
  m4 = m1;
  EXPECT_TRUE(m1 == m4);
}

TEST(RectangularBinaryMatrix, InitRaw) {
  const unsigned int nb_col = 80;
  uint64_t raw[nb_col];
  for(unsigned int i = 0; i < nb_col; ++i)
    raw[i] = random_bits();
  const RectangularBinaryMatrix m(raw, 19, nb_col);
  EXPECT_EQ(19u, m.r());
  EXPECT_EQ(80u, m.c());
  const uint64_t mask = ((uint64_t)1 << 19) - 1;
  for(unsigned int i = 0; i < nb_col; ++i)
    EXPECT_EQ(raw[i] & mask, m[i]);
}

TEST(RectangularBinaryMatrix, LowIdentity) {
  for(int r = 2; r < 64; r += 2) {
    for(int c = 2; c < 100; c += 2) {
      SCOPED_TRACE(::testing::Message() << "matrix " << r << "x" << c);
      RectangularBinaryMatrix m(r, c); // Matrix should be zeroed out
      mer_dna::k(c);
      mer_dna v;
      v.randomize();

      EXPECT_FALSE(m.is_low_identity());
      m.init_low_identity();
      EXPECT_TRUE(m.is_low_identity());

      uint64_t res = m.times(v);
      EXPECT_EQ(v.get_bits(0, std::min(r, c)), res);
    }
  }
}

/******************************
 * Matrix Vector product
 ******************************/
class MatrixVectorProd : public ::testing::TestWithParam< ::std::tr1::tuple<int, int> > {
public:
  unsigned int row, col;
  RectangularBinaryMatrix m;
  MatrixVectorProd() :
    row(::std::tr1::get<0>(GetParam())),
    col(::std::tr1::get<1>(GetParam())),
    m(row, col)
  {
    m.randomize(random_bits);
  }
};

TEST_P(MatrixVectorProd, Checks) {
  EXPECT_EQ(row, m.r());
  EXPECT_EQ(col, m.c());
}

TEST_P(MatrixVectorProd, AllOnes) {
  uint64_t v[2], res = 0;
  v[0] = v[1] = (uint64_t)-1;
  for(unsigned int i = 0; i < m.c(); ++i)
    res ^= m[i];
  EXPECT_EQ(res, m.times_loop(v));
#ifdef HAVE_INT128
  EXPECT_EQ(res, m.times_128(v));
#endif
#ifdef HAVE_SSE
  EXPECT_EQ(res, m.times_sse(v));
#endif
}

TEST_P(MatrixVectorProd, EveryOtherOnes) {
  uint64_t v[2], res = 0;
  v[0] = 0xaaaaaaaaaaaaaaaaUL;
  v[1] = 0xaaaaaaaaaaaaaaaaUL;
  for(unsigned int i = 0; i < m.c(); i += 2)
    res ^= m[i];
  EXPECT_EQ(res, m.times_loop(v));
#ifdef HAVE_INT128
  EXPECT_EQ(res, m.times_128(v));
#endif
#ifdef HAVE_SSE
  EXPECT_EQ(res, m.times_sse(v));
#endif

  v[0] >>= 1;
  v[1] >>= 1;
  res    = 0;
  for(unsigned int i = 1; i < m.c(); i += 2)
    res ^= m[i];
  EXPECT_EQ(res, m.times_loop(v));
#ifdef HAVE_INT128
  EXPECT_EQ(res, m.times_128(v));
#endif
#ifdef HAVE_SSE
  EXPECT_EQ(res, m.times_sse(v));
#endif
}

#if HAVE_SSE || HAVE_INT128
TEST_P(MatrixVectorProd, Optimizations) {
  static const int nb_tests = 100;
  const unsigned int nb_words = col / 64 + (col % 64 != 0);
  uint64_t v[nb_words];

  for(int i = 0; i < nb_tests; ++i) {
    // unsigned int r = 2 * (random() % 31 + 1);
    // unsigned int c = 2 * (random() % 100) + r;

    // RectangularBinaryMatrix m(r, c);
    // m.randomize(random_bits);

    for(unsigned int j = 0; j < nb_words; ++j)
      v[j] = random_bits();

    uint64_t res = m.times_loop((uint64_t*)v);
#ifdef HAVE_SSE
    EXPECT_EQ(res, m.times_sse((uint64_t*)v));
#endif
#ifdef HAVE_INT128
    EXPECT_EQ(res, m.times_128((uint64_t*)v));
#endif
  }
}
#endif // HAVE_SSE || HAVE_INT128

INSTANTIATE_TEST_CASE_P(MatrixVectorProdTest, MatrixVectorProd, ::testing::Combine(::testing::Range(1, 65, 1), // rows
                                                                                   ::testing::Range(2, 100, 2))); // cols

/******************************
 * Matrix product and inverse
 ******************************/
TEST(PseudoProduct, Dimensions) {
  RectangularBinaryMatrix m(30, 100), m1(32, 100), m2(30, 98);

  EXPECT_THROW(m.pseudo_multiplication(m1), std::domain_error);
  EXPECT_THROW(m.pseudo_multiplication(m2), std::domain_error);
}

TEST(PseudoProduct, Identity) {
  RectangularBinaryMatrix m(30, 100), i(30, 100);
  i.init_low_identity();
  m.randomize(random_bits);

  EXPECT_TRUE(i.pseudo_multiplication(m) == m);
}

TEST(PseudoProduct, Parity) {
  const unsigned int col_sizes[6] = { 50, 70, 126, 130, 64, 128 };
  const unsigned int nb_rows = 30;

  for(unsigned int k = 0; k < sizeof(col_sizes) / sizeof(unsigned int); ++k) {
    const unsigned int nb_cols = col_sizes[k];
    uint64_t *cols = new uint64_t[nb_cols];
    RectangularBinaryMatrix p(nb_rows, nb_cols);

    for(unsigned int j = 18; j < 19; ++j) {
      const uint64_t bits = ((uint64_t)1 << j) - 1;
      unsigned int i;
      for(i = 0; i < nb_cols; ++i)
        cols[i] = bits;
      RectangularBinaryMatrix m(cols, nb_rows, nb_cols);

      p = m.pseudo_multiplication(m);
      for(i = 0; i < nb_cols - nb_rows; ++i)
        EXPECT_EQ(__builtin_parity(bits) ? (uint64_t)0 : bits, p[i]);
      for( ; i < nb_cols; ++i)
        EXPECT_EQ(__builtin_parity(bits) ? bits : (uint64_t)0, p[i]);
    }
    delete [] cols;
  }
}

TEST(PseudoProduct, Inverse) {
  int full_rank = 0, singular = 0;
  for(unsigned int i = 0; i < 500; ++i) {
    unsigned int r = random() % 63 + 1;
    unsigned int c = 2 * ((random() % 100) + 1);
    SCOPED_TRACE(::testing::Message() << "Dimension " << r << "x" << c);
    RectangularBinaryMatrix m(r, c);
    m.randomize(random_bits);
    RectangularBinaryMatrix s(m);
    unsigned int rank = m.pseudo_rank();
    if(rank != c) {
      ++singular;
      EXPECT_THROW(m.pseudo_inverse(), std::domain_error);
    } else {
      ++full_rank;
      RectangularBinaryMatrix inv(m);
      EXPECT_NO_THROW(inv = m.pseudo_inverse());
      RectangularBinaryMatrix i = inv.pseudo_multiplication(m);
      EXPECT_TRUE(i.is_low_identity());
    }
    EXPECT_TRUE(s == m);
  }
  EXPECT_EQ(500, full_rank + singular);
  EXPECT_NE(0, full_rank);
}

TEST(PseudoProduct, Rank) {
  RectangularBinaryMatrix m(50, 100);
  for(unsigned int i = 0; i < 10; ++i) {
    m.randomize(random_bits);
    RectangularBinaryMatrix s(m);
    unsigned int rank = m.pseudo_rank();
    EXPECT_TRUE(rank <= m.c());
    EXPECT_TRUE(s == m);
  }
}

TEST(PseudoProduct, InitRandom) {
  RectangularBinaryMatrix m(50, 100);
  for(unsigned int i = 0; i < 10; ++i) {
    RectangularBinaryMatrix im(m.randomize_pseudo_inverse(random_bits));
    EXPECT_EQ((unsigned int)m.c(), m.pseudo_rank());
    EXPECT_EQ((unsigned int)m.c(), im.pseudo_rank());
    EXPECT_TRUE((m.pseudo_multiplication(im)).is_low_identity());
  }
}

TEST(PseudoProduct, VectorInverseMultiplication) {
  RectangularBinaryMatrix m(50, 132);
  RectangularBinaryMatrix im(m.randomize_pseudo_inverse());
  RectangularBinaryMatrix cim(im);
  EXPECT_TRUE(cim.pseudo_multiplication(m).is_low_identity());

  mer_dna::k(66);
  mer_dna v, iv, iv2;
  EXPECT_EQ(m.nb_words(), v.nb_words());

  for(int i = 0; i < 100; ++i) {
    SCOPED_TRACE(::testing::Message() << "i=" << i);
    v.randomize();
    uint64_t hash = m.times(v);
    iv = v; iv.set_bits(0, 50, hash);

    uint64_t lower = im.times(iv);
    EXPECT_EQ(v.get_bits(0, 50), lower);
  }
}


static const int speed_loop = 100000000;
TEST(MatrixProductSpeed, Loop) {
  RectangularBinaryMatrix m(50, 100);
  const unsigned int nb_words = m.c() / 64 + (m.c() % 64 != 0);
  uint64_t v[nb_words];
  for(unsigned int j = 0; j < nb_words; ++j)
    v[j] = random_bits();

  volatile uint64_t res = 0;
  for(int i = 0; i < speed_loop; ++i)
    res ^= m.times_loop((uint64_t*)v);
}

#ifdef HAVE_SSE
TEST(MatrixProductSpeed, SSE) {
  RectangularBinaryMatrix m(50, 100);
  const unsigned int nb_words = m.c() / 64 + (m.c() % 64 != 0);
  uint64_t v[nb_words];
  for(unsigned int j = 0; j < nb_words; ++j)
    v[j] = random_bits();

  volatile uint64_t res = 0;
  for(int i = 0; i < speed_loop; ++i)
    res ^= m.times_sse((uint64_t*)v);
}
#endif

#ifdef HAVE_INT128
TEST(MatrixProductSpeed, U128) {
  RectangularBinaryMatrix m(50, 100);
  const unsigned int nb_words = m.c() / 64 + (m.c() % 64 != 0);
  uint64_t v[nb_words];
  for(unsigned int j = 0; j < nb_words; ++j)
    v[j] = random_bits();

  volatile uint64_t res = 0;
  for(int i = 0; i < speed_loop; ++i)
    res ^= m.times_128((uint64_t*)v);
}
#endif

}
