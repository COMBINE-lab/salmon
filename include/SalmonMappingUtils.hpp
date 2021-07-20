/**
>HEADER
    Copyright (c) 2013-2019 Rob Patro rob@cs.umd.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#ifndef __SALMON_MAPPING_UTILS__
#define __SALMON_MAPPING_UTILS__

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iterator>
#include <vector>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"
#include "spdlog/spdlog.h"

// Future C++ convenience classes
#include "core/range.hpp"

#include "AlevinOpts.hpp"
#include "AlignmentGroup.hpp"
#include "ProgramOptionsGenerator.hpp"
#include "ReadExperiment.hpp"
#include "SalmonDefaults.hpp"
#include "SalmonMath.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"

#include "parallel_hashmap/phmap.h"
#include "pufferfish/MemChainer.hpp"
#include "pufferfish/MemCollector.hpp"
#include "pufferfish/PuffAligner.hpp"
#include "pufferfish/SAMWriter.hpp"
#include "pufferfish/SelectiveAlignmentUtils.hpp"
#include "pufferfish/Util.hpp"
#include "pufferfish/itlib/small_vector.hpp"
#include "pufferfish/ksw2pp/KSW2Aligner.hpp"
#include "pufferfish/metro/metrohash64.h"

namespace salmon {
namespace mapping_utils {

using MateStatus = pufferfish::util::MateStatus;
constexpr const int32_t invalid_score_ = std::numeric_limits<int32_t>::min();
constexpr const int32_t invalid_index_ = std::numeric_limits<int32_t>::min();
constexpr const size_t static_vec_size = 32;

class MappingScoreInfo {
public:
  MappingScoreInfo()
      : bestScore(invalid_score_), secondBestScore(invalid_score_),
        bestDecoyScore(invalid_score_), decoyThresh(1.0),
        collect_decoy_info_(false) {}

  MappingScoreInfo(double decoyThreshIn) : MappingScoreInfo() {
    decoyThresh = decoyThreshIn;
  }

  void collect_decoys(bool do_collect) { collect_decoy_info_ = do_collect; }
  bool collect_decoys() const { return collect_decoy_info_; }

  // clear everything but the decoy threshold
  void clear(size_t num_hits) {
    bestScore = invalid_score_;
    secondBestScore = invalid_score_;
    bestDecoyScore = invalid_score_;
    scores_.clear();
    scores_.resize(num_hits, invalid_score_);
    bestScorePerTranscript_.clear();
    perm_.clear();
    // NOTE: we do _not_ reset decoyThresh here
    if (collect_decoy_info_) {
      bool revert_to_static = best_decoy_hits.size() > static_vec_size;
      best_decoy_hits.clear();
      if (revert_to_static) {
        best_decoy_hits.revert_to_static();
      }
    }
  }

  bool haveOnlyDecoyMappings() const {
    // if the best non-decoy mapping has score less than decoyThresh *
    // bestDecoyScore and if the bestDecoyScore is a valid value, then we
    // have no valid non-decoy mappings.
    return (bestScore < static_cast<int32_t>(decoyThresh * bestDecoyScore)) and
           (bestDecoyScore >
            std::numeric_limits<decltype(bestDecoyScore)>::min());
  }

  inline bool update_decoy_mappings(int32_t hitScore, size_t idx,
                                    uint32_t tid) {
    const bool better_score = hitScore > bestDecoyScore;
    if (hitScore > bestDecoyScore) {
      bestDecoyScore = hitScore;
      if (collect_decoy_info_) {
        best_decoy_hits.clear();
        best_decoy_hits.push_back(std::make_pair(static_cast<int32_t>(idx),
                                                 static_cast<int32_t>(tid)));
      }
    } else if (collect_decoy_info_ and (hitScore == bestDecoyScore)) {
      best_decoy_hits.push_back(
          std::make_pair(static_cast<int32_t>(idx), static_cast<int32_t>(tid)));
    }
    return better_score;
  }

  int32_t bestScore;
  int32_t secondBestScore;
  int32_t bestDecoyScore;
  double decoyThresh;
  itlib::small_vector<std::pair<int32_t, int32_t>> best_decoy_hits;
  bool collect_decoy_info_;
  std::vector<int32_t> scores_;
  phmap::flat_hash_map<uint32_t, std::pair<int32_t, int32_t>>
      bestScorePerTranscript_;
  std::vector<std::pair<int32_t, int32_t>> perm_;
};

template <typename IndexT>
inline bool
initMapperSettings(SalmonOpts& salmonOpts, MemCollector<IndexT>& memCollector,
                   ksw2pp::KSW2Aligner& aligner,
                   pufferfish::util::AlignmentConfig& aconf,
                   pufferfish::util::MappingConstraintPolicy& mpol) {
  memCollector.configureMemClusterer(salmonOpts.maxOccsPerHit);
  double consensusFraction = (salmonOpts.consensusSlack == 0.0)
                                 ? 1.0
                                 : (1.0 - salmonOpts.consensusSlack);
  memCollector.setConsensusFraction(consensusFraction);
  memCollector.setHitFilterPolicy(salmonOpts.hitFilterPolicy);
  memCollector.setAltSkip(salmonOpts.mismatchSeedSkip);
  memCollector.setChainSubOptThresh(salmonOpts.pre_merge_chain_sub_thresh);

  // Initialize ksw aligner
  ksw2pp::KSW2Config config;
  // config.dropoff = -1;
  config.gapo = salmonOpts.gapOpenPenalty;
  config.gape = salmonOpts.gapExtendPenalty;
  config.bandwidth = salmonOpts.dpBandwidth;
  config.flag = 0;
  config.flag |= KSW_EZ_RIGHT;
  if (salmonOpts.computeCIGAR == false) {
        config.flag |= KSW_EZ_SCORE_ONLY;
  }
  int8_t a = static_cast<int8_t>(salmonOpts.matchScore);
  int8_t b = static_cast<int8_t>(salmonOpts.mismatchPenalty);
  ksw2pp::KSW2Aligner aligner2(static_cast<int8_t>(a), static_cast<int8_t>(b));
  aligner2.config() = config;
  std::swap(aligner, aligner2);

  aconf.refExtendLength = 20;
  aconf.fullAlignment = salmonOpts.fullLengthAlignment;
  aconf.mismatchPenalty = salmonOpts.mismatchPenalty;
  aconf.bestStrata = false;
  aconf.decoyPresent = false;
  aconf.matchScore = salmonOpts.matchScore;
  aconf.gapExtendPenalty = salmonOpts.gapExtendPenalty;
  aconf.gapOpenPenalty = salmonOpts.gapOpenPenalty;
  aconf.minScoreFraction = salmonOpts.minScoreFraction;
  aconf.mimicBT2 = salmonOpts.mimicBT2;
  aconf.mimicBT2Strict = salmonOpts.mimicStrictBT2;
  // aconf.allowOverhangSoftclip = salmonOpts.softclipOverhangs;
  aconf.allowSoftclip  = salmonOpts.softclip || salmonOpts.softclipOverhangs; 
  aconf.computeCIGAR = salmonOpts.computeCIGAR && (!salmonOpts.qmFileName.empty());
  aconf.endBonus = 5;
  aconf.end2end = !aconf.allowSoftclip;
  aconf.maxSoftclipFraction = salmonOpts.maxSoftclipFraction;
  aconf.useAlignmentCache = !salmonOpts.disableAlignmentCache;
  // aconf.alignmentMode = pufferfish::util::PuffAlignmentMode::SCORE_ONLY;

  // we actually care about the softclips in the cigar string 
  // if we are writing output and softclipping (or softclipping of overhangs) is enabled
  // if ( (!salmonOpts.qmFileName.empty()) and (salmonOpts.softclip or salmonOpts.softclipOverhangs) ) {
  //   aconf.alignmentMode = pufferfish::util::PuffAlignmentMode::APPROXIMATE_CIGAR;
  // }

  mpol.noOrphans = !salmonOpts.allowOrphans;
  // TODO : PF_INTEGRATION
  // decide how we want to set this
  // I think we don't want to allow general discordant reads
  // e.g. both map to same strand or map too far away, but
  // we want the "allowDovetail" option to determine if
  // a dovetail read is considered concordant or discordant
  mpol.noDiscordant = true;
  mpol.noDovetail = !salmonOpts.allowDovetail;
  aconf.noDovetail = mpol.noDovetail;

  mpol.setPostMergeChainSubThresh(salmonOpts.post_merge_chain_sub_thresh);
  mpol.setOrphanChainSubThresh(salmonOpts.orphan_chain_sub_thresh);

  return true;
}

inline void updateRefMappings(uint32_t tid, int32_t hitScore, bool isCompat,
                              size_t idx,
                              const std::vector<Transcript>& transcripts,
                              int32_t invalidScore,
                              salmon::mapping_utils::MappingScoreInfo& msi) {
  auto& scores = msi.scores_;
  scores[idx] = hitScore;
  auto& t = transcripts[tid];
  bool isDecoy = t.isDecoy();
  double decoyCutoff =
      static_cast<int32_t>(msi.decoyThresh * msi.bestDecoyScore);

  // if (hitScore < decoyCutoff or (hitScore == invalidScore)) { }

  if (isDecoy) {
    // NOTE: decide here if we need to process any of this if the
    // current score is < the best (non-decoy) score. I think not.

    // if this is a decoy and its score is better than the best decoy score
    bool did_update = msi.update_decoy_mappings(hitScore, idx, tid);
    (void)did_update;
    return;
  } else if (hitScore < decoyCutoff or (hitScore == invalidScore)) {
    // if the current score is to a valid target but doesn't
    // exceed the necessary decoy threshold, then skip it.
    return;
  }
  // otherwise, we have a "high-scoring" hit to a non-decoy

  auto& perm = msi.perm_;
  auto& bestScorePerTranscript = msi.bestScorePerTranscript_;
  // removing duplicate hits from a read to the same transcript
  auto it = bestScorePerTranscript.find(tid);
  if (it == bestScorePerTranscript.end()) {
    // if we didn't have any alignment for this transcript yet, then
    // this is the current best
    bestScorePerTranscript[tid].first = hitScore;
    bestScorePerTranscript[tid].second = idx;
  } else if ((hitScore > it->second.first) or
             (hitScore == it->second.first and isCompat)) {
    // otherwise, if we had an alignment for this transcript and it's
    // better than the current best, then set the best score to this
    // alignment's score, and invalidate the previous alignment
    it->second.first = hitScore;
    scores[it->second.second] = invalidScore;
    it->second.second = idx;
  } else {
    // otherwise, there is already a better mapping for this transcript.
    scores[idx] = invalidScore;
  }

  if (hitScore > msi.bestScore) {
    msi.secondBestScore = msi.bestScore;
    msi.bestScore = hitScore;
  }
  perm.push_back(std::make_pair(idx, tid));
}

inline void filterAndCollectAlignments(
    std::vector<pufferfish::util::JointMems>& jointHits, uint32_t readLen,
    uint32_t mateLen, bool singleEnd, bool tryAlign, bool hardFilter,
    double scoreExp, double minAlnProb,
    salmon::mapping_utils::MappingScoreInfo& msi,
    std::vector<pufferfish::util::QuasiAlignment>& jointAlignments) {

  auto invalidScore = std::numeric_limits<decltype(msi.bestDecoyScore)>::min();
  if (msi.bestDecoyScore == invalidScore) {
    msi.bestDecoyScore = invalidScore + 1;
  }

  int32_t decoyThreshold = static_cast<decltype(msi.bestDecoyScore)>(
      msi.decoyThresh * msi.bestDecoyScore);

  // auto filterScore = (bestDecoyScore < secondBestScore) ? secondBestScore :
  // bestDecoyScore;
  auto& scores = msi.scores_;
  auto& perm = msi.perm_;

  // throw away any pairs for which we should not produce valid alignments :
  // ======
  // If we are doing soft-filtering (default), we remove those not exceeding the
  // bestDecoyScore If we are doing hard-filtering, we remove those less than
  // the bestScore
  perm.erase(std::remove_if(
                 perm.begin(), perm.end(),
                 [&scores, hardFilter, &msi, decoyThreshold, invalidScore](
                     const std::pair<int32_t, int32_t>& idxtid) -> bool {
                   return !hardFilter ? scores[idxtid.first] < decoyThreshold
                                      : scores[idxtid.first] < msi.bestScore;
                   // return !hardFilter ?  scores[idxtid.first] < filterScore :
                   // scores[idxtid.first] < bestScore;
                 }),
             perm.end());

  // Unlike RapMap, pufferfish doesn't guarantee the hits computed above are in
  // order by transcript, so we find the permutation of indices that puts things
  // in transcript order.
  std::sort(perm.begin(), perm.end(),
            [](const std::pair<int32_t, int32_t>& p1,
               const std::pair<int32_t, int32_t>& p2) {
              return p1.second < p2.second;
            });

  // moving our alinged / score jointMEMs over to QuasiAlignment objects
  double bestScoreD = static_cast<double>(msi.bestScore);
  for (auto& idxTxp : perm) {
    int32_t ctr = idxTxp.first;
    int32_t tid = idxTxp.second;
    auto& jointHit = jointHits[ctr];

    double currScore = scores[ctr];
    double v = bestScoreD - currScore;
    // why -1?
    double estAlnProb = hardFilter ? -1.0 : std::exp(-scoreExp * v);
    // skip any alignment with aln prob < minAlnProb
    if (!hardFilter and (estAlnProb < minAlnProb)) {
      continue;
    }

    if (singleEnd or jointHit.isOrphan()) {
      readLen = jointHit.isLeftAvailable() ? readLen : mateLen;
      jointAlignments.emplace_back(
          tid,                                        // reference id
          jointHit.orphanClust()->getTrFirstHitPos(), // reference pos
          jointHit.orphanClust()->isFw,               // fwd direction
          readLen,                                    // read length
          jointHit.orphanClust()->cigar,              // cigar string
          jointHit.fragmentLen,                       // fragment length
          false);
      auto& qaln = jointAlignments.back();
      // NOTE : score should not be filled in from a double
      qaln.score = !tryAlign
                       ? static_cast<int32_t>(jointHit.orphanClust()->coverage)
                       : jointHit.alignmentScore;
      qaln.estAlnProb(estAlnProb);
      // NOTE : wth is numHits?
      qaln.numHits =
          static_cast<uint32_t>(jointHits.size()); // orphanClust()->coverage;
      qaln.mateStatus = jointHit.mateStatus;
      if (singleEnd) {
        qaln.mateLen = readLen;
        qaln.mateCigar.clear();
        qaln.matePos = 0;
        qaln.mateIsFwd = true;
        qaln.mateScore = 0;
        qaln.mateStatus = MateStatus::SINGLE_END;
      }
    } else {
      jointAlignments.emplace_back(
          tid,                                    // reference id
          jointHit.leftClust->getTrFirstHitPos(), // reference pos
          jointHit.leftClust->isFw,               // fwd direction
          readLen,                                // read length
          jointHit.leftClust->cigar,              // cigar string
          jointHit.fragmentLen,                   // fragment length
          true);                                  // properly paired
      // Fill in the mate info
      auto& qaln = jointAlignments.back();
      qaln.mateLen = mateLen;
      qaln.mateCigar = jointHit.rightClust->cigar;
      qaln.matePos =
          static_cast<int32_t>(jointHit.rightClust->getTrFirstHitPos());
      qaln.mateIsFwd = jointHit.rightClust->isFw;
      qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
      // NOTE : wth is numHits?
      qaln.numHits = static_cast<uint32_t>(jointHits.size());
      // NOTE : score should not be filled in from a double
      qaln.score = !tryAlign
                       ? static_cast<int32_t>(jointHit.leftClust->coverage)
                       : jointHit.alignmentScore;
      qaln.estAlnProb(estAlnProb);
      qaln.mateScore = !tryAlign
                           ? static_cast<int32_t>(jointHit.rightClust->coverage)
                           : jointHit.mateAlignmentScore;
    }
  }
  // done moving our alinged / score jointMEMs over to QuasiAlignment objects
}

inline void filterAndCollectAlignmentsDecoy(
    std::vector<pufferfish::util::JointMems>& jointHits, uint32_t readLen,
    uint32_t mateLen, bool singleEnd, bool tryAlign, bool hardFilter,
    double scoreExp, double minAlnProb,
    salmon::mapping_utils::MappingScoreInfo& msi,
    std::vector<pufferfish::util::QuasiAlignment>& jointAlignments) {
  // NOTE: this function should only be called in the case that we have valid
  // decoy mappings to report. Currently, this happens only when there are *no
  // valid non-decoy* mappings. Further, this function will only add equally
  // *best* decoy mappings to the output jointAlignments object regardless of
  // the the status of hardFilter (i.e. no sub-optimal decoy mappings will be
  // reported).
  (void)hardFilter;
  (void)minAlnProb;
  (void)scoreExp;
  double estAlnProb = 1.0; // std::exp(-scoreExp * 0.0);
  for (auto& idxTxp : msi.best_decoy_hits) {
    int32_t ctr = idxTxp.first;
    int32_t tid = idxTxp.second;
    auto& jointHit = jointHits[ctr];

    if (singleEnd or jointHit.isOrphan()) {
      readLen = jointHit.isLeftAvailable() ? readLen : mateLen;
      jointAlignments.emplace_back(
          tid,                                        // reference id
          jointHit.orphanClust()->getTrFirstHitPos(), // reference pos
          jointHit.orphanClust()->isFw,               // fwd direction
          readLen,                                    // read length
          jointHit.orphanClust()->cigar,              // cigar string
          jointHit.fragmentLen,                       // fragment length
          false);
      auto& qaln = jointAlignments.back();
      // NOTE : score should not be filled in from a double
      qaln.score = !tryAlign
                       ? static_cast<int32_t>(jointHit.orphanClust()->coverage)
                       : jointHit.alignmentScore;
      qaln.estAlnProb(estAlnProb);
      // NOTE : wth is numHits?
      qaln.numHits =
          static_cast<uint32_t>(jointHits.size()); // orphanClust()->coverage;
      qaln.mateStatus = jointHit.mateStatus;
      if (singleEnd) {
        qaln.mateLen = readLen;
        qaln.mateCigar.clear();
        qaln.matePos = 0;
        qaln.mateIsFwd = true;
        qaln.mateScore = 0;
        qaln.mateStatus = MateStatus::SINGLE_END;
      }
    } else {
      jointAlignments.emplace_back(
          tid,                                    // reference id
          jointHit.leftClust->getTrFirstHitPos(), // reference pos
          jointHit.leftClust->isFw,               // fwd direction
          readLen,                                // read length
          jointHit.leftClust->cigar,              // cigar string
          jointHit.fragmentLen,                   // fragment length
          true);                                  // properly paired
      // Fill in the mate info
      auto& qaln = jointAlignments.back();
      qaln.mateLen = mateLen;
      qaln.mateCigar = jointHit.rightClust->cigar;
      qaln.matePos =
          static_cast<int32_t>(jointHit.rightClust->getTrFirstHitPos());
      qaln.mateIsFwd = jointHit.rightClust->isFw;
      qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
      // NOTE : wth is numHits?
      qaln.numHits = static_cast<uint32_t>(jointHits.size());
      // NOTE : score should not be filled in from a double
      qaln.score = !tryAlign
                       ? static_cast<int32_t>(jointHit.leftClust->coverage)
                       : jointHit.alignmentScore;
      qaln.estAlnProb(estAlnProb);
      qaln.mateScore = !tryAlign
                           ? static_cast<int32_t>(jointHit.rightClust->coverage)
                           : jointHit.mateAlignmentScore;
    }
  } // end for over best decoy hits
}

namespace pasc {

constexpr int32_t invalid_frag_len = std::numeric_limits<int32_t>::min();
constexpr int32_t invalid_mate_pos = std::numeric_limits<int32_t>::min();

/* from pesc
struct simple_hit {
    bool is_fw{false};
    bool mate_is_fw{false};
    int32_t pos{-1};
    float score{0.0};
    uint32_t num_hits{0};
    uint32_t tid{std::numeric_limits<uint32_t>::max()};
    int32_t mate_pos{std::numeric_limits<int32_t>::max()};
    int32_t fragment_length{std::numeric_limits<int32_t>::max()};
    inline bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over)
{ int32_t signed_txp_len = static_cast<int32_t>(txp_len); return (pos >
-max_over) and ((pos + read_len) < (signed_txp_len + max_over));
    }
    inline bool has_mate() const { return mate_pos != invalid_mate_pos; }
    inline bool mate_is_mapped() const { return mate_pos != invalid_mate_pos; }
    inline int32_t frag_len() const {
        return (fragment_length != invalid_frag_len) ? fragment_length : 0;
    }
};
*/
struct simple_hit {
  bool is_fw{false};
  int32_t pos{-1};
  float score{0.0};
  uint32_t num_hits{0};
  uint32_t tid{std::numeric_limits<uint32_t>::max()};
  // int32_t mate_pos{std::numeric_limits<int32_t>::max()};
  // int32_t fragment_length{std::numeric_limits<int32_t>::max()};

  // UNPAIRED_LEFT, UNPAIRED_RIGHT, PAIRED_FR, PAIRED_RF
  PairingStatus pairing_status{PairingStatus::UNPAIRED_RIGHT};

  bool valid_pos(int32_t read_len, uint32_t txp_len, int32_t max_over) {
    int32_t signed_txp_len = static_cast<int32_t>(txp_len);
    return (pos > -max_over) and
           ((pos + read_len) < (signed_txp_len + max_over));
  }

  bool operator<(const simple_hit& hit2) {
    if (tid != hit2.tid) { // hit1 and hit2 are on different transcripts
      return tid < hit2.tid;
    } else { // hit1 and hit2 are on the same transcript
      // NOTE: Is this what we actually want?
      if (is_fw) { // fw < rc
        return true;
      } else {
        return false;
      }
    }
  }
};

enum class HitDirection : uint8_t { FW, RC, BOTH };

/*
struct sketch_hit_info {
  // add a hit to the current target that occurs in the forward
  // orientation with respect to the target.
  bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k,
              int32_t max_stretch, float score_inc) {
    (void)rl;
    (void)k;
    bool added{false};

    // since hits are collected by moving _forward_ in the
    // read, if this is a fw hit, it should be moving
    // forward in the reference. Only add it if this is
    // the case.  This ensure that we don't
    // double-count a k-mer that might occur twice on
    // this target.
    if (ref_pos > last_ref_pos_fw and read_pos > last_read_pos_fw) {
      if (last_read_pos_fw == -1) {
        approx_pos_fw = ref_pos - read_pos;
      } else {
        if ((ref_pos - approx_pos_fw) > max_stretch) {
          return false;
        }
      }
      // if (last_ref_pos_fw > -1 and (ref_pos > last_ref_pos_fw + 15)) { return
      // false; }
      last_ref_pos_fw = ref_pos;
      last_read_pos_fw = read_pos;
      fw_score += score_inc;
      ++fw_hits;
      added = true;
    }
    return added;
  }

  // add a hit to the current target that occurs in the forward
  // orientation with respect to the target.
  bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k,
              int32_t max_stretch, float score_inc) {

    bool added{false};
    // since hits are collected by moving _forward_ in the
    // read, if this is an rc hit, it should be moving
    // backwards in the reference. Only add it if this is
    // the case.
    // This ensures that we don't double-count a k-mer that
    // might occur twice on this target.
    if (ref_pos < last_ref_pos_rc and read_pos > last_read_pos_rc) {
      approx_pos_rc = (ref_pos - (rl - (read_pos + k)));
      if (last_read_pos_rc == -1) {
        approx_end_pos_rc = ref_pos + read_pos;
      } else {
        if (approx_end_pos_rc - approx_pos_rc > max_stretch) {
          return false;
        }
      }
      // if (last_ref_pos_rc > -1 and ref_pos < last_ref_pos_rc - 15) { return
      // false; }
      last_ref_pos_rc = ref_pos;
      last_read_pos_rc = read_pos;
      rc_score += score_inc;
      ++rc_hits;
      added = true;
    }
    return added;
  }

  inline uint32_t max_hits_for_target() { return std::max(fw_hits, rc_hits); }

  // true if forward, false if rc
  // second element is score
  inline HitDirection best_hit_direction() {
    int32_t fw_minus_rc =
        static_cast<int32_t>(fw_hits) - static_cast<int32_t>(rc_hits);
    return (fw_minus_rc > 0)
               ? HitDirection::FW
               : ((fw_minus_rc < 0) ? HitDirection::RC : HitDirection::BOTH);
  }

  inline simple_hit get_fw_hit(PairingStatus ps) {
    return simple_hit{true,
                      approx_pos_fw,
                      fw_score,
                      fw_hits,
                      std::numeric_limits<uint32_t>::max(),
                      ps};
  }

  inline simple_hit get_rc_hit(PairingStatus ps) {
    return simple_hit{false,
                      approx_pos_rc,
                      rc_score,
                      rc_hits,
                      std::numeric_limits<uint32_t>::max(),
                      ps};
  }

  inline std::string to_string() {
    std::stringstream ss;
    ss << "fw_hits: " << fw_hits << ", fw_score : " << fw_score
       << ", fw_pos : " << approx_pos_fw << " || rc_hits: " << rc_hits
       << ", rc_score: " << rc_score << ", rc_pos: " << approx_pos_rc;
    return ss.str();
  }

  int32_t last_read_pos_fw{-1};
  int32_t last_read_pos_rc{-1};

  int32_t last_ref_pos_fw{-1};
  int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};

  int32_t approx_pos_fw{-1};
  int32_t approx_pos_rc{-1};
  int32_t approx_end_pos_rc{-1};

  uint32_t fw_hits{0};
  uint32_t rc_hits{0};
  float fw_score{0.0};
  float rc_score{0.0};
};*/

struct sketch_hit_info {
    // add a hit to the current target that occurs in the forward
    // orientation with respect to the target.
    bool add_fw(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch,
                float score_inc) {
        (void)rl;
        (void)k;
        bool added{false};

        // since hits are collected by moving _forward_ in the
        // read, if this is a fw hit, it should be moving
        // forward in the reference. Only add it if this is
        // the case.  This ensure that we don't
        // double-count a k-mer that might occur twice on
        // this target.
        if (ref_pos > last_ref_pos_fw and read_pos > last_read_pos_fw) {
            if (last_read_pos_fw == -1) {
                approx_pos_fw = ref_pos - read_pos;
            } else {
                if ((ref_pos - approx_pos_fw) > max_stretch) { return false; }
            }
            last_ref_pos_fw = ref_pos;
            last_read_pos_fw = read_pos;
            fw_score += score_inc;
            ++fw_hits;
            added = true;
        }
        return added;
    }

    // add a hit to the current target that occurs in the forward
    // orientation with respect to the target.
    bool add_rc(int32_t ref_pos, int32_t read_pos, int32_t rl, int32_t k, int32_t max_stretch,
                float score_inc) {
        bool added{false};
        // since hits are collected by moving _forward_ in the
        // read, if this is an rc hit, it should be moving
        // backwards in the reference. 
        // In general, we only add the hit if this is the case.
        // This ensures that we don't double-count a k-mer that
        // might occur twice on this target.

        // we have a special case here; what if the same exact 
        // k-mer (i.e. not just the same sequence but same position
        // on the query) occurs more than one time on this refernece?
        //
        // In that case, the GENERAL case code will already have 
        // processed and seen a k-mer with the read position 
        // equal to `read_pos`.  In the case below, we see 
        // a hit with the *same* read pos again (one or more times).
        // 
        // Here, we swap out the previous hit having this read_pos 
        // if the position of the current hit on the read is
        // the same and the position on the reference is greater 
        // (this is a heuristic to help in the case of tandem repeats or 
        // highly-repetitive subsequence).
        // NOTE: consider if a similar heuristic should be
        // adopted for the forward case.
        if ((read_pos == last_read_pos_rc) and (ref_pos > last_ref_pos_rc) and
            (ref_pos < rightmost_bound_rc)) {

            last_ref_pos_rc = ref_pos;
            // if the read_pos was the same as the first read pos
            // then also update the approx_end_pos_rc accordingly
            // NOTE: for the time being don't mess with this position
            // empirically this does better, but if we really want 
            // to optimize this for accuracy we need a better general
            // heuristic.
            // if (read_pos == first_read_pos_rc) {
            //   approx_end_pos_rc = ref_pos + read_pos;
            //}

            return added;
        }

        // GENERAL case
        if (ref_pos < last_ref_pos_rc and read_pos > last_read_pos_rc) {
            approx_pos_rc = (ref_pos - (rl - (read_pos + k)));
            if (last_read_pos_rc == -1) {
                approx_end_pos_rc = ref_pos + read_pos;
                first_read_pos_rc = read_pos;
            } else {
                if (approx_end_pos_rc - approx_pos_rc > max_stretch) { return false; }
            }
            rc_score += score_inc;
            ++rc_hits;

            // new
            rightmost_bound_rc = last_ref_pos_rc;

            last_ref_pos_rc = ref_pos;
            last_read_pos_rc = read_pos;
            added = true;
        }
        return added;
    }

    inline uint32_t max_hits_for_target() { return std::max(fw_hits, rc_hits); }

    // true if forward, false if rc
    // second element is score
    inline HitDirection best_hit_direction() {
        int32_t fw_minus_rc = static_cast<int32_t>(fw_hits) - static_cast<int32_t>(rc_hits);
        return (fw_minus_rc > 0) ? HitDirection::FW
                                 : ((fw_minus_rc < 0) ? HitDirection::RC : HitDirection::BOTH);
    }

    inline simple_hit get_fw_hit(PairingStatus ps) {
        return simple_hit{true, approx_pos_fw,
                          fw_score, fw_hits, std::numeric_limits<uint32_t>::max(), ps};
    }

    inline simple_hit get_rc_hit(PairingStatus ps) {
        return simple_hit{false, approx_pos_rc,
                          rc_score, rc_hits, std::numeric_limits<uint32_t>::max(), ps};
    }

    inline std::string to_string() {
        std::stringstream ss;
        ss << "fw_hits: " << fw_hits << ", fw_score : " << fw_score
           << ", fw_pos : " << approx_pos_fw << " || rc_hits: " << rc_hits
           << ", rc_score: " << rc_score << ", rc_pos: " << approx_pos_rc;
        return ss.str();
    }

    int32_t last_read_pos_fw{-1};
    int32_t last_read_pos_rc{-1};
    int32_t rightmost_bound_rc{std::numeric_limits<int32_t>::max()};

    // marks the read position (key) of the 
    // first hit we see in the rc direction
    int32_t first_read_pos_rc{-1};

    int32_t last_ref_pos_fw{-1};
    int32_t last_ref_pos_rc{std::numeric_limits<int32_t>::max()};

    int32_t approx_pos_fw{-1};
    int32_t approx_pos_rc{-1};
    int32_t approx_end_pos_rc{-1};

    uint32_t fw_hits{0};
    uint32_t rc_hits{0};
    float fw_score{0.0};
    float rc_score{0.0};
};

/**
 * This contains the information necessary to collect hits and perform
 * mapping on a read using PASC.
 **/
template <typename IndexT> struct mapping_cache_info {
public:
  mapping_cache_info(IndexT* midx, size_t max_occ_default_in,
                     size_t max_occ_recover_in, size_t max_read_occ_in)
      : idx(midx), k(midx->k()), hit_searcher(midx),
        max_occ_default(max_occ_default_in),
        max_occ_recover(max_occ_recover_in),
        attempt_occ_recover(max_occ_recover > max_occ_default),
        max_read_occ(max_read_occ_in), alt_max_occ(max_read_occ_in) {}

  inline void clear() {
    map_type = salmon::utils::MappingType::UNMAPPED;
    hit_searcher.clear();
    hit_map.clear();
    accepted_hits.clear();
    alt_max_occ = max_read_occ;
    has_matching_kmers = false;
  }

  // will store how the read mapped
  salmon::utils::MappingType map_type{salmon::utils::MappingType::UNMAPPED};

  // map from reference id to hit info
  phmap::flat_hash_map<uint32_t, sketch_hit_info> hit_map;

  // where the mappings that pass will be stored
  std::vector<simple_hit> accepted_hits;

  // map to recall the number of unmapped reads we see
  // for each barcode
  phmap::flat_hash_map<uint64_t, uint32_t> unmapped_bc_map;

  // pointer to the underlying index
  IndexT* idx{nullptr};

  // the k-mer size of our index
  size_t k{0};

  // implements the PASC algorithm
  MemCollector<IndexT> hit_searcher;

  // the number of occurrences of a hit above which we
  // won't consider it
  size_t max_occ_default;

  // the number of occurences of a hit above which we
  // won't consider it, even in recovery mode
  size_t max_occ_recover;

  // will we attempt recovery if there are no hits with
  // frequency below max_occ_default?
  bool attempt_occ_recover;

  // maximum number of places a read can occur and still
  // be considered.
  size_t max_read_occ;

  size_t alt_max_occ;

  // to perform queries
  pufferfish::util::QueryCache qc;

  // regardless of having full mappings, did any k-mers match
  bool has_matching_kmers{false};
};

/**
 * This function will map the read given by `read_seq`
 * using PASC. The relevant parameters (e.g. maximum ambiguity
 * of seed and maximum number of allowable mappings) are stored 
 * in the `map_cache` structure.  The `PairingStatus` is passed in 
 * by the caller and designates if the read being mapped is the 
 * left end or the right end.
 **/
template <typename IndexT>
inline bool map_read(std::string* read_seq,
                     mapping_cache_info<IndexT>& map_cache, PairingStatus ps,
                     bool verbose = false) {
  map_cache.clear();
  // rebind map_cache variables to
  // local names
  IndexT* qidx = map_cache.idx;
  auto& qc = map_cache.qc;
  auto& hit_searcher = map_cache.hit_searcher;
  auto& hit_map = map_cache.hit_map;
  auto& accepted_hits = map_cache.accepted_hits;
  auto& map_type = map_cache.map_type;
  const bool attempt_occ_recover = map_cache.attempt_occ_recover;
  auto k = map_cache.k;
  int32_t signed_k = static_cast<int32_t>(k);

  // collect the set of matching seeds
  map_cache.has_matching_kmers =
      hit_searcher.get_raw_hits_sketch(*read_seq, qc, true, false);
  bool early_stop = false;

  // if there were hits
  if (map_cache.has_matching_kmers) {
    uint32_t num_valid_hits{0};
    uint64_t total_occs{0};
    uint64_t largest_occ{0};
    auto& raw_hits = hit_searcher.get_left_hits();

    // SANITY
    decltype(raw_hits[0].first) prev_read_pos = -1;
    // the maximum span the supporting k-mers of a
    // mapping position are allowed to have.
    // NOTE this is still > read_length b/c the stretch is measured wrt the
    // START of the terminal k-mer.
    int32_t max_stretch = static_cast<int32_t>(read_seq->length() * 1.0);

    // a raw hit is a pair of read_pos and a projected hit

    // the least frequent hit for this fragment.
    uint64_t min_occ = std::numeric_limits<uint64_t>::max();

    // this is false by default and will be set to true
    // if *every* collected hit for this fragment occurs
    // max_occ_default times or more.
    bool had_alt_max_occ = false;
    int32_t signed_rl = static_cast<int32_t>(read_seq->length());
    auto collect_mappings_from_hits =
        [&max_stretch, &min_occ, &hit_map, &num_valid_hits, &total_occs,
         &largest_occ, &early_stop, signed_rl, k, signed_k, qidx,
         verbose](auto& raw_hits, auto& prev_read_pos, auto& max_allowed_occ,
                  auto& had_alt_max_occ) -> bool {
      for (auto& raw_hit : raw_hits) {
        auto& read_pos = raw_hit.first;
        auto& proj_hits = raw_hit.second;
        auto& refs = proj_hits.refRange;

        // Keep track of the *least* ambiguous hit we see.
        // If all hits occur more than `max_allowed_occ` times,
        // but the least ambiguous hit (which occurs `min_occ`
        // times) has < `max_occ_recover` occurrences, then 
        // we will use all hits occurring `min_occ` number 
        // of times or less to try top recover the read's mappings.
        uint64_t num_occ = static_cast<uint64_t>(refs.size());
        min_occ = std::min(min_occ, num_occ);
        had_alt_max_occ = true;

        // we visit every hit that occurs the allowable number of time
        bool still_have_valid_target = false;
        prev_read_pos = read_pos;
        if (num_occ <= max_allowed_occ) {
          total_occs += num_occ;
          largest_occ = (num_occ > largest_occ) ? num_occ : largest_occ;
          float score_inc = 1.0;

          // vist each mapped reference position (mrp) where this hit 
          // occurs.
          for (auto& pos_it : refs) {
            const auto& ref_pos_ori = proj_hits.decodeHit(pos_it);
            uint32_t tid =
                static_cast<uint32_t>(qidx->getRefId(pos_it.transcript_id()));
            int32_t pos = static_cast<int32_t>(ref_pos_ori.pos);
            bool ori = ref_pos_ori.isFW;
            auto& target = hit_map[tid];

            // Why >= here instead of == ?
            // Because hits can happen on the same target in both the forward
            // and rc orientations, it is possible that we start the loop with
            // the target having num_valid_hits hits in a given orientation (o)
            // we see a new hit for this target in oriention o (now it has
            // num_valid_hits + 1) then we see a hit for this target in
            // orientation rc(o).  We still want to add / consider this hit, but
            // max_hits_for_target() > num_valid_hits. So, we must allow for
            // that here.
            if (target.max_hits_for_target() >= num_valid_hits) {
              if (ori) {
                target.add_fw(pos, static_cast<int32_t>(read_pos), signed_rl,
                              signed_k, max_stretch, score_inc);
              } else {
                target.add_rc(pos, static_cast<int32_t>(read_pos), signed_rl,
                              signed_k, max_stretch, score_inc);
              }

              // if this target is still valid after evaluating this hit 
              // then we will have at least one globally valid target. If 
              // no targets are valid after evaluating all mrps for this 
              // hit then we can exit the search early.
              still_have_valid_target |=
                  (target.max_hits_for_target() >= num_valid_hits + 1);
            }

          } // DONE: for (auto &pos_it : refs)

          ++num_valid_hits;

          // if there are no targets reaching the valid hit threshold, then
          // break early
          if (!still_have_valid_target) {
            return true;
          }

        } // DONE : if (num_occ <= max_allowed_occ)
      }   // DONE : for (auto& raw_hit : raw_hits)

      return false;
    };

    bool _discard = false;
    auto mao_first_pass = map_cache.max_occ_default - 1;
    early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos,
                                            mao_first_pass, _discard);

    // If our default threshold was too stringent, then fallback to a more
    // liberal threshold and look up the k-mers that occur the least frequently.
    // Specifically, if the min occuring hits have frequency < max_occ_recover
    // (2500 by default) times, then collect the min occuring hits to get the
    // mapping.
    if (attempt_occ_recover and (min_occ >= map_cache.max_occ_default) and
        (min_occ < map_cache.max_occ_recover)) {
      prev_read_pos = -1;
      uint64_t max_allowed_occ = min_occ;
      early_stop = collect_mappings_from_hits(raw_hits, prev_read_pos,
                                              max_allowed_occ, had_alt_max_occ);
    }

    uint32_t best_alt_hits = 0;

    for (auto& kv : hit_map) {
      auto best_hit_dir = kv.second.best_hit_direction();
      // if the best direction is FW or BOTH, add the fw hit
      // otherwise add the RC.
      auto mapping = (best_hit_dir != HitDirection::RC)
                         ? kv.second.get_fw_hit(ps)
                         : kv.second.get_rc_hit(ps);

      if (mapping.num_hits >= num_valid_hits) {
        mapping.tid = kv.first;
        accepted_hits.emplace_back(mapping);
        // if we had equally good hits in both directions
        // add the rc hit here (since we added the fw)
        // above if the best hit was either FW or BOTH
        if (best_hit_dir == HitDirection::BOTH) {
          auto second_hit = kv.second.get_rc_hit(ps);
          second_hit.tid = kv.first;
          accepted_hits.emplace_back(second_hit);
        }
      } else {
        // best_alt_score = mapping.score > best_alt_score ? mapping.score :
        // best_alt_score;
        best_alt_hits =
            mapping.num_hits > best_alt_hits ? mapping.num_hits : best_alt_hits;
      }
    }

    map_cache.alt_max_occ =
        had_alt_max_occ ? accepted_hits.size() : map_cache.max_occ_default;

    /*
     * This rule; if enabled, allows through mappings missing a single hit, if
    there
     * was no mapping with all hits. NOTE: this won't work with the current
    early-exit
     * optimization however.
     if (accepted_hits.empty() and (num_valid_hits > 1) and (best_alt_hits >=
     num_valid_hits
     - 1)) { for (auto& kv : hit_map) { auto mapping = kv.second.get_best_hit();
    if (mapping.num_hits >= best_alt_hits) {
    //if (mapping.valid_pos(signed_read_len, transcripts[kv.first].RefLength,
    10)) { mapping.tid = kv.first; accepted_hits.emplace_back(mapping);
    //}
    }
    }
    }
    */
  } // DONE : if (rh)

  // If the read mapped to > maxReadOccs places, discard it
  if (accepted_hits.size() > map_cache.alt_max_occ) {
    accepted_hits.clear();
    map_type = salmon::utils::MappingType::UNMAPPED;
  } else if (!accepted_hits.empty()) {
    map_type = salmon::utils::MappingType::SINGLE_MAPPED;
  }

  return early_stop;
}

} // namespace pasc

} // namespace mapping_utils
} // namespace salmon

#endif //__SALMON_MAPPING_UTILS__
