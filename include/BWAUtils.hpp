#ifndef __BWA_UTILS_HPP__
#define __BWA_UTILS_HPP__

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
}

namespace bwautils {

// Function modified from bwt_smem1a:
// https://github.com/lh3/bwa/blob/eb428d7d31ced059ad39af2701a22ebe6d175657/bwt.c#L289
/**
 * Search for the k-mer of length @len starting at @q.
 * Return true if an interval is found for the k-mer and false
 * otherwise. The appropriate bwt interval will be placed
 * in @resInterval upon success.
 *
 */
bool getIntervalForKmer(const bwt_t* bwt, // the bwt index
                        int len,          // k-mer length
                        const uint8_t* q, // query
                        bwtintv_t& resInterval);

// NOTE: $max_intv is not currently used in BWA-MEM
// NOTE: Modified from the original functions to take an initial interval for
// the search query
int bwt_smem1a_with_kmer(const bwt_t* bwt, int len, const uint8_t* q, int x,
                         int min_intv, uint64_t max_intv,
                         bwtintv_t initial_interval, bwtintv_v* mem,
                         bwtintv_v* tmpvec[2]);

int bwt_smem1_with_kmer(const bwt_t* bwt, int len, const uint8_t* q, int x,
                        int min_intv, bwtintv_t initial_interval,
                        bwtintv_v* mem, bwtintv_v* tmpvec[2]);
} // namespace bwautils

#endif // __BWA_UTILS_HPP__
