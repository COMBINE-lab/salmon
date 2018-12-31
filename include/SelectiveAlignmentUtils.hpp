#ifndef __SELECTIVE_ALIGNMENT_UTILS__
#define __SELECTIVE_ALIGNMENT_UTILS__


#include "ksw2pp/KSW2Aligner.hpp"
#include "metro/metrohash64.h"
#include "tsl/hopscotch_map.h"

namespace selective_alignment {
  namespace utils {

    /// Get alignment score
    namespace salmon {
      namespace mapping {
        struct CacheEntry {
          uint64_t hashKey;
          int32_t alnScore;
          CacheEntry(uint64_t hashKeyIn, int32_t alnScoreIn) :
            hashKey(hashKeyIn), alnScore(alnScoreIn) {}
        };
      }
    }

    // Use a passthrough hash for the alignment cache, because
    // the key *is* the hash.
    namespace salmon {
      namespace hashing {
        struct PassthroughHash {
          std::size_t operator()(uint64_t const& u) const { return u; }
        };
      }
    }

    using AlnCacheMap = tsl::hopscotch_map<uint64_t, int32_t, salmon::hashing::PassthroughHash>;
    //tsl::robin_map<uint64_t, int32_t, salmon::hashing::PassthroughHash>;
    //std::unordered_map<uint64_t, int32_t, salmon::hashing::PassthroughHash>;

inline int32_t getAlnScore(
                           ksw2pp::KSW2Aligner& aligner,
                           ksw_extz_t& ez,
                           int32_t pos, const char* rptr, int32_t rlen,
                           char* tseq, int32_t tlen,
                           int8_t mscore,
                           int8_t mmcost,
                           int32_t maxScore,
                           rapmap::utils::ChainStatus chainStat,
                           bool multiMapping, // was there > 1 hit for this read
                           bool mimicStrictBT2,
                           uint32_t buf,
                           AlnCacheMap& alnCache) {
  // If this was a perfect match, don't bother to align or compute the score
  if (chainStat == rapmap::utils::ChainStatus::PERFECT) {
    return maxScore;
  }

  auto ungappedAln = [mscore, mmcost](char* ref, const char* query, int32_t len) -> int32_t {
    int32_t ungappedScore{0};
    for (int32_t i = 0; i < len; ++i) {
      char c1 = *(ref + i);
      char c2 = *(query + i);
      c1 = (c1 == 'N' or c2 == 'N') ? c2 : c1;
      ungappedScore += (c1 == c2) ? mscore : mmcost;
    }
    return ungappedScore;
  };

  int32_t s{std::numeric_limits<int32_t>::lowest()};
  // TODO : Determine what is the most "appropriate" penalty for
  // an overhang (based on the scoring function).
  bool invalidStart = (pos < 0);
  if (invalidStart) { rptr += -pos; rlen += pos; pos = 0; }

  // if we are trying to mimic Bowtie2 with RSEM params
  if (invalidStart and mimicStrictBT2) { return s; }

  if (pos < tlen) {
    bool doUngapped{(!invalidStart) and (chainStat == rapmap::utils::ChainStatus::UNGAPPED)};
    buf = (doUngapped) ? 0 : buf;
    auto lnobuf = static_cast<uint32_t>(tlen - pos);
    auto lbuf = static_cast<uint32_t>(rlen+buf);
    auto useBuf = (lbuf < lnobuf);
    uint32_t tlen1 = std::min(lbuf, lnobuf);
    char* tseq1 = tseq + pos;
    ez.max_q = ez.max_t = ez.mqe_t = ez.mte_q = -1;
    ez.max = 0, ez.mqe = ez.mte = KSW_NEG_INF;
    ez.n_cigar = 0;

    uint64_t hashKey{0};
    bool didHash{false};
    if (!alnCache.empty()) {
      // hash the reference string
      uint32_t keyLen = useBuf ? tlen1 - buf : tlen1;
      MetroHash64::Hash(reinterpret_cast<uint8_t*>(tseq1), keyLen, reinterpret_cast<uint8_t*>(&hashKey), 0);
      didHash = true;
      // see if we have this hash
      auto hit = alnCache.find(hashKey);
      // if so, we know the alignment score
      if (hit != alnCache.end()) {
        s = hit->second;
      }
    }
    // If we got here with s == -1, we don't have the score cached
    if (s == std::numeric_limits<int32_t>::lowest()) {
      if (doUngapped) {
        // signed version of tlen1
        int32_t tlen1s = tlen1;
        int32_t alnLen = rlen < tlen1s ? rlen : tlen1s;
        s = ungappedAln(tseq1, rptr, alnLen);
      } else {
        aligner(rptr, rlen, tseq1, tlen1, &ez, ksw2pp::EnumToType<ksw2pp::KSW2AlignmentType::EXTENSION>());
        s = std::max(ez.mqe, ez.mte);
      }

      if (multiMapping) { // don't bother to fill up a cache unless this is a multi-mapping read
        if (!didHash) {
          uint32_t keyLen = useBuf ? tlen1 - buf : tlen1;
          MetroHash64::Hash(reinterpret_cast<uint8_t*>(tseq1), keyLen, reinterpret_cast<uint8_t*>(&hashKey), 0);
        }
        alnCache[hashKey] = s;
      } // was a multi-mapper

    }
  }
  return s;
}

  } // namespace utils
} // namespace selective_alignment

#endif //__SELECTIVE_ALIGNMENT_UTILS__
