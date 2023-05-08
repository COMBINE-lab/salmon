#include "ksw2pp/KSW2Aligner.hpp"
#include <limits>
#include <iostream>

/*
extern void ksw_extz2_sse2(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
extern void ksw_extz2_sse41(void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat, int8_t q, int8_t e, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez);
*/

namespace ksw2pp {

// The following is from ksw2_dispatch.c : https://github.com/lh3/minimap2/blob/master/ksw2_dispatch.c
#define SIMD_SSE     0x1
#define SIMD_SSE2    0x2
#define SIMD_SSE3    0x4
#define SIMD_SSSE3   0x8
#define SIMD_SSE4_1  0x10
#define SIMD_SSE4_2  0x20
#define SIMD_AVX     0x40
#define SIMD_AVX2    0x80
#define SIMD_AVX512F 0x100

#ifdef KSW_USE_ARM
#else
#ifndef _MSC_VER
// adapted from https://github.com/01org/linux-sgx/blob/master/common/inc/internal/linux/cpuid_gnu.h
void __cpuidex(int cpuid[4], int func_id, int subfunc_id)
{
#if defined(__x86_64__)
	asm volatile ("cpuid"
			: "=a" (cpuid[0]), "=b" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#else // on 32bit, ebx can NOT be used as PIC code
	asm volatile ("xchgl %%ebx, %1; cpuid; xchgl %%ebx, %1"
			: "=a" (cpuid[0]), "=r" (cpuid[1]), "=c" (cpuid[2]), "=d" (cpuid[3])
			: "0" (func_id), "2" (subfunc_id));
#endif
}
#endif

int x86_simd(void)
{
	int flag = 0, cpuid[4], max_id;
	__cpuidex(cpuid, 0, 0);
	max_id = cpuid[0];
	if (max_id == 0) return 0;
	__cpuidex(cpuid, 1, 0);
	if (cpuid[3]>>25&1) flag |= SIMD_SSE;
	if (cpuid[3]>>26&1) flag |= SIMD_SSE2;
	if (cpuid[2]>>0 &1) flag |= SIMD_SSE3;
	if (cpuid[2]>>9 &1) flag |= SIMD_SSSE3;
	if (cpuid[2]>>19&1) flag |= SIMD_SSE4_1;
	if (cpuid[2]>>20&1) flag |= SIMD_SSE4_2;
	if (cpuid[2]>>28&1) flag |= SIMD_AVX;
	if (max_id >= 7) {
		__cpuidex(cpuid, 7, 0);
		if (cpuid[1]>>5 &1) flag |= SIMD_AVX2;
		if (cpuid[1]>>16&1) flag |= SIMD_AVX512F;
	}
	return flag;
}
#endif // KSW_USE_ARM
// end of ksw2_dispatch.c here

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
    #ifdef KSW_USE_ARM
    haveSSE41 = true;
    haveSSE2 = true;
    #else
    unsigned int simd = x86_simd();
    haveSSE41 = (simd & SIMD_SSE4_1);
    haveSSE2 = (simd & SIMD_SSE2);
    #endif // KSW_USE_ARM
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
    for (j = 0; j < m; ++j) {
      mat_[(m - 1) * m + j] = 0;
    }
}

KSW2Aligner::KSW2Aligner(std::vector<int8_t> mat) {
    #ifdef KSW_USE_ARM
    haveSSE41 = true;
    haveSSE2 = true;
    #else
    unsigned int simd = x86_simd();
    haveSSE41 = (simd & SIMD_SSE4_1);
    haveSSE2 = (simd & SIMD_SSE2);
    #endif // KSW_USE_ARM

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
int
KSW2Aligner::transformSequenceKSW2(const char* const queryOriginal, const int queryLength,
                                   std::vector<unsigned char>& queryTransformed) {
  // Alphabet is constructed from letters that are present in sequences.
  // Each letter is assigned an ordinal number, starting from 0 up to
  // alphabetLength - 1, and new query and target are created in which letters
  // are replaced with their ordinal numbers. This query and target are used in
  // all the calculations later.
  queryTransformed.resize(queryLength, 0);

  int i = 0;
  for (int i = 0; i < queryLength; i++) {
    uint8_t c = static_cast<uint8_t>(queryOriginal[i]);
    queryTransformed[i] = seq_nt4_table_loc[c];
  }
  return 4;
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
int
KSW2Aligner::transformSequencesKSW2(const char* const queryOriginal, const int queryLength,
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
                            int cutoff,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  // NOTE: all ksw extension aligner calls clear out ez 
  // *internally*.  This is why we do not need to (and do not)
  // clear it out here.
  // ksw_reset_extz(ez);
  auto qlen = queryLength;
  auto tlen = targetLength;
  int asize = transformSequencesKSW2(queryOriginal, queryLength, targetOriginal,
                                    targetLength, query_, target_);
  (void)asize;
  int8_t q = config_.gapo;
  int8_t e = config_.gape;
  int max_qt_len = (queryLength > targetLength) ? queryLength : targetLength;
  int w = (config_.bandwidth > max_qt_len) ? max_qt_len : config_.bandwidth;
  int z = config_.dropoff;
  if (haveSSE41) {
    ksw_extz2_sse41(kalloc_allocator_.get(), qlen, query_.data(), tlen,
                target_.data(), config_.alphabetSize, mat_.data(), q, e, w, z,
                config_.end_bonus, config_.flag, ez, cutoff);
  } else if (haveSSE2) {
    ksw_extz2_sse2(kalloc_allocator_.get(), qlen, query_.data(), tlen,
                  target_.data(), config_.alphabetSize, mat_.data(), q, e, w, z,
                  config_.end_bonus, config_.flag, ez, cutoff);
  } else {
    std::abort();
  }
  return ez->score;
}

int KSW2Aligner::operator()(const char* const queryOriginal,
                            const int queryLength,
                            const char* const targetOriginal,
                            const int targetLength,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  return this->operator()(queryOriginal, queryLength, targetOriginal,
                          targetLength, &result_, std::numeric_limits<int>::min(),
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
                           targetLength, &result_, std::numeric_limits<int>::min(),
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
  ksw_reset_extz(ez);
  auto qlen = queryLength;
  auto tlen = targetLength;
  int asize = transformSequencesKSW2(queryOriginal, queryLength, targetOriginal,
                                    targetLength, query_, target_);
  (void)asize;
  int q = config_.gapo;
  int e = config_.gape;
  int max_qt_len = (queryLength > targetLength) ? queryLength : targetLength;
  //int w = (config_.bandwidth > max_qt_len) ? max_qt_len : config_.bandwidth;
  int w = max_qt_len;

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
  ksw_reset_extz(ez);
  auto qlen = queryLength;
  auto tlen = targetLength;
  int q = config_.gapo;
  int e = config_.gape;
  int max_qt_len = (queryLength > targetLength) ? queryLength : targetLength;
  //int w = (config_.bandwidth > max_qt_len) ? max_qt_len : config_.bandwidth;
  int w = max_qt_len;

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
    ret = this->operator()(query_, queryLength, target_, targetLength, &result_, std::numeric_limits<int>::min(),
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
                            const int targetLength, ksw_extz_t* ez, int cutoff,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  // NOTE: all ksw extension aligner calls clear out ez 
  // *internally*.  This is why we do not need to (and do not)
  // clear it out here.
  //ksw_reset_extz(ez);
  auto qlen = queryLength;
  auto tlen = targetLength;
  int8_t q = config_.gapo;
  int8_t e = config_.gape;
  int max_qt_len = (queryLength > targetLength) ? queryLength : targetLength;
  int w = (config_.bandwidth > max_qt_len) ? max_qt_len : config_.bandwidth;
  int z = config_.dropoff;
  if (haveSSE41) {
    ksw_extz2_sse41(kalloc_allocator_.get(), qlen, query_, tlen, target_,
                config_.alphabetSize, mat_.data(), q, e, w, z, config_.end_bonus, config_.flag,
                ez, cutoff);
  } else if (haveSSE2) {
    ksw_extz2_sse2(kalloc_allocator_.get(), qlen, query_, tlen, target_,
                  config_.alphabetSize, mat_.data(), q, e, w, z, config_.end_bonus, config_.flag,
                  ez, cutoff);
  } else {
    std::abort();
  }
  return ez->score;
}

int KSW2Aligner::operator()(const uint8_t* const query_, const int queryLength,
                            const uint8_t* const target_,
                            const int targetLength,
                            EnumToType<KSW2AlignmentType::EXTENSION>) {
  return this->operator()(query_, queryLength, target_, targetLength, &result_, std::numeric_limits<int>::min(),
                          EnumToType<KSW2AlignmentType::EXTENSION>());
}

} // namespace ksw2pp
