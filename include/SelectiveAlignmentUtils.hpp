#ifndef __SELECTIVE_ALIGNMENT_UTILS__
#define __SELECTIVE_ALIGNMENT_UTILS__


#include "ksw2pp/KSW2Aligner.hpp"
#include "metro/metrohash64.h"
#include "tsl/hopscotch_map.h"
#include "edlib.h"

namespace selective_alignment {
  namespace utils {

    enum class AlignmentPolicy : uint8_t { DEFAULT, BT2, BT2_STRICT };

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


    inline bool recoverOrphans(std::string& leftRead,
                               std::string& rightRead,
                               std::string& rc1,
                               std::string& rc2,
                               const std::vector<Transcript>& transcripts,
                               const std::vector<rapmap::utils::QuasiAlignment>& leftHits,
                               const std::vector<rapmap::utils::QuasiAlignment>& rightHits,
                               std::vector<rapmap::utils::QuasiAlignment>& jointHits) {
      using QuasiAlignment = rapmap::utils::QuasiAlignment;
      auto* r1 = leftRead.data();
      auto* r2 = rightRead.data();
      auto l1 = static_cast<int32_t>(leftRead.length());
      auto l2 = static_cast<int32_t>(rightRead.length());
      // We compute the reverse complements below only if we
      // need them and don't have them.
      char* r1rc = nullptr;
      char* r2rc = nullptr;

      const char* windowSeq = nullptr;
      int32_t windowLength = -1;

      int32_t maxDistRight = l2 / 4;
      int32_t maxDistLeft = l1 / 4;
      constexpr const int32_t signedZero{0};
      int32_t lreadLen = l1;
      int32_t rreadLen = l2;

      for (auto && lh : leftHits) {
        auto txpID = lh.tid;
        auto& txp = transcripts[lh.tid];
        const char* tseq = txp.Sequence();

        int32_t startPos = -1;
        const char* rptr = r2;
        auto rlen = l2;
        auto maxDist = maxDistRight;

        auto leftChainStatus = lh.chainStatus.getLeft();
        auto rightChainStatus = rapmap::utils::ChainStatus::REGULAR;

        int32_t lpos{lh.allPositions.front()};
        int32_t rpos{-1};

        bool lfwd = lh.fwd;
        bool rfwd = !lfwd;

        // if this hit is forward, look downstream, else upstream
        if ( lh.fwd ) {
          if (!r2rc){
            rapmap::utils::reverseRead(rightRead, rc2);
            r2rc = const_cast<char*>(rc2.data());
          }
          rptr = r2rc;
          startPos = std::max(signedZero, static_cast<int32_t>(lpos));
          windowLength = std::min(1000, static_cast<int32_t>(txp.RefLength - startPos));
        } else {
          rptr = r2;
          int32_t endPos = std::min(static_cast<int32_t>(txp.RefLength), static_cast<int32_t>(lpos) + l1);
          startPos = std::max(signedZero,  endPos - 1000);
          windowLength = std::min(1000, endPos);
        }
        windowSeq = tseq + startPos;

        EdlibAlignResult result = edlibAlign(rptr, rlen, windowSeq, windowLength,
                                             edlibNewAlignConfig(maxDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
        if (result.editDistance > -1) {
          rpos = startPos + result.endLocations[0] - l2;
          // If we consider only a single position per transcript
          int32_t startRead1 = std::max(lpos, signedZero);
          int32_t startRead2 = std::max(rpos, signedZero);
          bool read1First{(startRead1 < startRead2)};
          int32_t fragStartPos = read1First ? startRead1 : startRead2;
          int32_t fragEndPos = read1First ?
            (startRead2 + l2) : (startRead1 + l1);
          uint32_t fragLen = fragEndPos - fragStartPos;
          jointHits.emplace_back(txpID,
                                 lpos,
                                 lfwd,
                                 lreadLen,
                                 fragLen, true);
          // Fill in the mate info
          auto& qaln = jointHits.back();
          qaln.mateLen = rlen;
          qaln.matePos = rpos;
          qaln.mateIsFwd = rfwd;
          jointHits.back().mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
          jointHits.back().chainStatus = rapmap::utils::FragmentChainStatus(leftChainStatus, rightChainStatus);
        }
        edlibFreeAlignResult(result);
      }

      for (auto && rh : rightHits) {
        auto txpID = rh.tid;
        auto& txp = transcripts[rh.tid];
        const char* tseq = txp.Sequence();

        int32_t startPos = -1;
        const char* rptr = r1;
        auto rlen = l1;
        auto maxDist = maxDistLeft;

        auto leftChainStatus = rapmap::utils::ChainStatus::REGULAR;
        auto rightChainStatus = rh.chainStatus.getRight();

        int32_t lpos{-1};
        int32_t rpos{rh.allPositions.front()};

        bool lfwd = !rh.fwd;
        bool rfwd = rh.fwd;

        // if this hit is forward, look downstream, else upstream
        if ( rh.fwd ) {
          if (!r1rc){
            rapmap::utils::reverseRead(leftRead, rc1);
            r1rc = const_cast<char*>(rc1.data());
          }
          rptr = r1rc;
          startPos = std::max(signedZero, static_cast<int32_t>(rpos));
          windowLength = std::min(1000, static_cast<int32_t>(txp.RefLength - startPos));
        } else {
          rptr = r1;
          int32_t endPos = std::min(static_cast<int32_t>(txp.RefLength), static_cast<int32_t>(rpos) + l2);
          startPos = std::max(signedZero,  endPos - 1000);
          windowLength = std::min(1000, endPos);
        }
        windowSeq = tseq + startPos;

        EdlibAlignResult result = edlibAlign(rptr, rlen, windowSeq, windowLength,
                                             edlibNewAlignConfig(maxDist, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE));
        if (result.editDistance > -1) {
          lpos = startPos + result.endLocations[0] - l1;
          // If we consider only a single position per transcript
          int32_t startRead1 = std::max(lpos, signedZero);
          int32_t startRead2 = std::max(rpos, signedZero);
          bool read1First{(startRead1 < startRead2)};
          int32_t fragStartPos = read1First ? startRead1 : startRead2;
          int32_t fragEndPos = read1First ?
            (startRead2 + l2) : (startRead1 + l1);
          uint32_t fragLen = fragEndPos - fragStartPos;
          jointHits.emplace_back(txpID,
                                 lpos,
                                 lfwd,
                                 lreadLen,
                                 fragLen, true);
          // Fill in the mate info
          auto& qaln = jointHits.back();
          qaln.mateLen = rlen;
          qaln.matePos = rpos;
          qaln.mateIsFwd = rfwd;
          jointHits.back().mateStatus = rapmap::utils::MateStatus::PAIRED_END_PAIRED;
          jointHits.back().chainStatus = rapmap::utils::FragmentChainStatus(leftChainStatus, rightChainStatus);
        }
        edlibFreeAlignResult(result);
      }

      std::sort(jointHits.begin(), jointHits.end(), [](const QuasiAlignment& a, const QuasiAlignment& b) {
          return a.tid < b.tid;
        });


      return true;
    }


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
                           AlignmentPolicy ap,
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
  bool invalidEnd = (pos + rlen >= tlen);
  if (invalidStart) { rptr += -pos; rlen += pos; pos = 0; }

  // if we are trying to mimic Bowtie2 with RSEM params
  if (invalidStart or invalidEnd) {
    switch (ap) {
    case AlignmentPolicy::BT2:
    case AlignmentPolicy::BT2_STRICT:
      return s;
    case AlignmentPolicy::DEFAULT:
    default:
      break;
    }
  }

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
        /*
        auto startBuf = std::min(pos, static_cast<int32_t>(buf));
        int32_t bpos = pos - startBuf;
        char* tseqTemp = tseq + bpos;
        uint32_t tlenTemp = tlen1 + startBuf;
        EdlibAlignResult result = edlibAlign(rptr, rlen, tseqTemp, tlenTemp,
                                             edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_LOC));
        auto spos = result.startLocations[0];
        tseq1 = tseq + bpos + spos;
        */
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
