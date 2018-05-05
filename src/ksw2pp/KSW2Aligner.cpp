#include "KSW2Aligner.hpp"
#include "ksw2.h"

namespace ksw2pp {
unsigned char seq_nt4_table_loc[256] = {
    0, 1, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 4, 1, 4, 4, 4, 2,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};

KSW2Aligner::KSW2Aligner(int8_t match, int8_t mismatch) {
  query_.clear();
  target_.clear();
  kalloc_allocator_.reset(km_init());
  int a = match;
  int b = mismatch;
  int m = 5;
  int i, j;
  // use the ``simple'' scoring matrix for now.
  mat_.resize(m * m);

  a = a < 0 ? -a : a;
  b = b > 0 ? -b : b;
  for (i = 0; i < m - 1; ++i) {
    for (j = 0; j < m - 1; ++j)
      mat_[i * m + j] = i == j ? a : b;
    mat_[i * m + m - 1] = 0;
  }
  for (j = 0; j < m; ++j)
    mat_[(m - 1) * m + j] = 0;
}

KSW2Aligner::KSW2Aligner(std::vector<int8_t> mat) {
  query_.clear();
  target_.clear();
  kalloc_allocator_.reset(km_init());
  mat_ = mat;
}

/**
 * from https://github.com/rob-p/edlib/blob/read-aligner/edlib/src/edlib.cpp
 * Takes char query and char target, recognizes alphabet and transforms them
 * into unsigned char sequences where elements in sequences are not any more
 * letters of alphabet, but their index in alphabet. Most of internal edlib
 * functions expect such transformed sequences. This function will allocate
 * queryTransformed and targetTransformed, so make sure to free them when done.
 * Example:
 *   Original sequences: "ACT" and "CGT".
 *   Alphabet would be recognized as ['A', 'C', 'T', 'G']. Alphabet length = 4.
 *   Transformed sequences: [0, 1, 2] and [1, 3, 2].
 * @param [in] queryOriginal
 * @param [in] queryLength
 * @param [in] targetOriginal
 * @param [in] targetLength
 * @param [out] queryTransformed  It will contain values in range [0, alphabet
 * length - 1].
 * @param [out] targetTransformed  It will contain values in range [0, alphabet
 * length - 1].
 * @return  Alphabet length - number of letters in recognized alphabet.
 */
static int
transformSequencesDNA(const char* const queryOriginal, const int queryLength,
                      const char* const targetOriginal, const int targetLength,
                      std::vector<unsigned char>& queryTransformed,
                      std::vector<unsigned char>& targetTransformed) {
  // Alphabet is constructed from letters that are present in sequences.
  // Each letter is assigned an ordinal number, starting from 0 up to
  // alphabetLength - 1, and new query and target are created in which letters
  // are replaced with their ordinal numbers. This query and target are used in
  // all the calculations later.
  queryTransformed.resize(queryLength, 0);
  targetTransformed.resize(targetLength, 0);

  int i = 0;
  int lengthDiff = targetLength - queryLength;
  int m = queryLength < targetLength ? queryLength : targetLength;
  int M = lengthDiff >= 0 ? (m + lengthDiff) : (m - lengthDiff);
  for (; i < m; i++) {
    uint8_t c = static_cast<uint8_t>(queryOriginal[i]);
    uint8_t c2 = static_cast<uint8_t>(targetOriginal[i]);
    queryTransformed[i] = seq_nt4_table_loc[c];
    targetTransformed[i] = seq_nt4_table_loc[c2];
  }
  if (lengthDiff > 0) {
    for (; i < M; i++) {
      uint8_t c = static_cast<uint8_t>(targetOriginal[i]);
      targetTransformed[i] = seq_nt4_table_loc[c];
    }
  } else if (lengthDiff < 0) {
    for (; i < M; i++) {
      uint8_t c = static_cast<uint8_t>(queryOriginal[i]);
      queryTransformed[i] = seq_nt4_table_loc[c];
    }
  }
  return 4;
}

int KSW2Aligner::operator()(const char* const queryOriginal,
                            const int queryLength,
                            const char* const targetOriginal,
                            const int targetLength, ksw_extz_t* ez,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  // auto ez = &result_;
  auto qlen = queryLength;
  auto tlen = targetLength;
  int asize = transformSequencesDNA(queryOriginal, queryLength, targetOriginal,
                                    targetLength, query_, target_);
  (void)asize;
  int q = config_.gapo;
  int e = config_.gape;
  int w = config_.bandwidth;
  int z = config_.dropoff;
  ksw_extz2_sse(kalloc_allocator_.get(), qlen, query_.data(), tlen,
                target_.data(), config_.alphabetSize, mat_.data(), q, e, w, z,
                config_.flag, ez);
  return ez->score;
}

int KSW2Aligner::operator()(const char* const queryOriginal,
                            const int queryLength,
                            const char* const targetOriginal,
                            const int targetLength,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  return this->operator()(queryOriginal, queryLength, targetOriginal,
                          targetLength, &result_,
                          EnumToType<KSW2AlignmentType::EXTENSION>());
}

int KSW2Aligner::operator()(const char* const queryOriginal,
                            const int queryLength,
                            const char* const targetOriginal,
                            const int targetLength) {
  int ret{0};
  switch (config_.atype) {
  case KSW2AlignmentType::EXTENSION:
    ret = this->operator()(queryOriginal, queryLength, targetOriginal,
                           targetLength, &result_,
                           EnumToType<KSW2AlignmentType::EXTENSION>());
    break;
  case KSW2AlignmentType::GLOBAL:
    ret = this->operator()(queryOriginal, queryLength, targetOriginal,
                           targetLength, &result_,
                           EnumToType<KSW2AlignmentType::GLOBAL>());
    break;
  }
  return ret;
}

int KSW2Aligner::operator()(const char* const queryOriginal,
                            const int queryLength,
                            const char* const targetOriginal,
                            const int targetLength, ksw_extz_t* ez,
                            EnumToType<KSW2AlignmentType::GLOBAL>) {
  // auto ez = &result_;
  auto qlen = queryLength;
  auto tlen = targetLength;
  int asize = transformSequencesDNA(queryOriginal, queryLength, targetOriginal,
                                    targetLength, query_, target_);
  (void)asize;
  int q = config_.gapo;
  int e = config_.gape;
  int w = config_.bandwidth;
  ez->score =
      (config_.flag & KSW_EZ_SCORE_ONLY)
          ? ksw_gg2(kalloc_allocator_.get(), qlen, query_.data(), tlen,
                    target_.data(), config_.alphabetSize, mat_.data(), q, e, w,
                    0, 0, 0)
          : ksw_gg2_sse(kalloc_allocator_.get(), qlen, query_.data(), tlen,
                        target_.data(), config_.alphabetSize, mat_.data(), q, e,
                        w, &ez->m_cigar, &ez->n_cigar, &ez->cigar);
  return ez->score;
}

int KSW2Aligner::operator()(const char* const queryOriginal,
                            const int queryLength,
                            const char* const targetOriginal,
                            const int targetLength,
                            EnumToType<KSW2AlignmentType::GLOBAL>) {
  return this->operator()(queryOriginal, queryLength, targetOriginal,
                          targetLength, &result_,
                          EnumToType<KSW2AlignmentType::GLOBAL>());
}

int KSW2Aligner::operator()(const uint8_t* const query_, const int queryLength,
                            const uint8_t* const target_,
                            const int targetLength, ksw_extz_t* ez,
                            EnumToType<KSW2AlignmentType::GLOBAL>) {
  auto qlen = queryLength;
  auto tlen = targetLength;
  int q = config_.gapo;
  int e = config_.gape;
  int w = config_.bandwidth;
  ez->score =
      (config_.flag & KSW_EZ_SCORE_ONLY)
          ? ksw_gg2(kalloc_allocator_.get(), qlen, query_, tlen, target_,
                    config_.alphabetSize, mat_.data(), q, e, w, 0, 0, 0)
          : ksw_gg2_sse(kalloc_allocator_.get(), qlen, query_, tlen, target_,
                        config_.alphabetSize, mat_.data(), q, e, w,
                        &ez->m_cigar, &ez->n_cigar, &ez->cigar);
  return ez->score;
}

int KSW2Aligner::operator()(const uint8_t* const query_, const int queryLength,
                            const uint8_t* const target_,
                            const int targetLength,
                            EnumToType<KSW2AlignmentType::GLOBAL>) {
  return this->operator()(query_, queryLength, target_, targetLength, &result_,
                          EnumToType<KSW2AlignmentType::GLOBAL>());
}

int KSW2Aligner::operator()(const uint8_t* const query_, const int queryLength,
                            const uint8_t* const target_,
                            const int targetLength) {
  int ret{0};
  switch (config_.atype) {
  case KSW2AlignmentType::EXTENSION:
    ret = this->operator()(query_, queryLength, target_, targetLength, &result_,
                           EnumToType<KSW2AlignmentType::EXTENSION>());
    break;
  case KSW2AlignmentType::GLOBAL:
    ret = this->operator()(query_, queryLength, target_, targetLength, &result_,
                           EnumToType<KSW2AlignmentType::GLOBAL>());
    break;
  }
  return ret;
}

int KSW2Aligner::operator()(const uint8_t* const query_, const int queryLength,
                            const uint8_t* const target_,
                            const int targetLength, ksw_extz_t* ez,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  auto qlen = queryLength;
  auto tlen = targetLength;
  int q = config_.gapo;
  int e = config_.gape;
  int w = config_.bandwidth;
  int z = config_.dropoff;
  ksw_extz2_sse(kalloc_allocator_.get(), qlen, query_, tlen, target_,
                config_.alphabetSize, mat_.data(), q, e, w, z, config_.flag,
                ez);
  return ez->score;
}

int KSW2Aligner::operator()(const uint8_t* const query_, const int queryLength,
                            const uint8_t* const target_,
                            const int targetLength,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  return this->operator()(query_, queryLength, target_, targetLength, &result_,
                          EnumToType<KSW2AlignmentType::EXTENSION>());
}

} // namespace ksw2pp
