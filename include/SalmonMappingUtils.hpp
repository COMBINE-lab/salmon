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

#include "SalmonDefaults.hpp"
#include "SalmonMath.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"
#include "AlignmentGroup.hpp"
#include "ProgramOptionsGenerator.hpp"
#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"

#include "pufferfish/Util.hpp"
#include "pufferfish/MemCollector.hpp"
#include "pufferfish/MemChainer.hpp"
#include "pufferfish/SAMWriter.hpp"
#include "pufferfish/PuffAligner.hpp"
#include "pufferfish/ksw2pp/KSW2Aligner.hpp"
#include "pufferfish/metro/metrohash64.h"
#include "pufferfish/SelectiveAlignmentUtils.hpp"
#include "pufferfish/chobo/small_vector.hpp"
#include "parallel_hashmap/phmap.h"

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
            bestDecoyScore(invalid_score_), decoyThresh(1.0), collect_decoy_info_(false) {}

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
        scores_.clear(); scores_.resize(num_hits, invalid_score_);
        bestScorePerTranscript_.clear();
        perm_.clear();
        // NOTE: we do _not_ reset decoyThresh here
        if (collect_decoy_info_) {
          bool revert_to_static = best_decoy_hits.size() > static_vec_size;
          best_decoy_hits.clear();
          if (revert_to_static) { best_decoy_hits.revert_to_static(); }
        }
      }

      bool haveOnlyDecoyMappings() const {
        // if the best non-decoy mapping has score less than decoyThresh *
        // bestDecoyScore and if the bestDecoyScore is a valid value, then we
        // have no valid non-decoy mappings.
        return (bestScore < static_cast<int32_t>(decoyThresh * bestDecoyScore)) and
               (bestDecoyScore > std::numeric_limits<decltype(bestDecoyScore)>::min());
      }
      
      inline bool update_decoy_mappings(int32_t hitScore, size_t idx, uint32_t tid) {
        const bool better_score = hitScore > bestDecoyScore;
        if (hitScore > bestDecoyScore) {
          bestDecoyScore = hitScore;
          if (collect_decoy_info_) { 
            best_decoy_hits.clear();
            best_decoy_hits.push_back(std::make_pair(static_cast<int32_t>(idx), static_cast<int32_t>(tid)));
          }
        } else if (collect_decoy_info_ and (hitScore == bestDecoyScore)){
          best_decoy_hits.push_back(std::make_pair(static_cast<int32_t>(idx), static_cast<int32_t>(tid)));
        }
        return better_score;
      }

      int32_t bestScore;
      int32_t secondBestScore;
      int32_t bestDecoyScore;
      double decoyThresh;
      chobo::small_vector<std::pair<int32_t, int32_t>> best_decoy_hits;
      bool collect_decoy_info_;
      std::vector<int32_t> scores_;
      phmap::flat_hash_map<uint32_t, std::pair<int32_t, int32_t>> bestScorePerTranscript_;
      std::vector<std::pair<int32_t, int32_t>> perm_;
    };

template <typename IndexT>
inline bool initMapperSettings(SalmonOpts& salmonOpts, MemCollector<IndexT>& memCollector, ksw2pp::KSW2Aligner& aligner,
                        pufferfish::util::AlignmentConfig& aconf, pufferfish::util::MappingConstraintPolicy& mpol) {
  memCollector.configureMemClusterer(salmonOpts.maxOccsPerHit);
  double consensusFraction = (salmonOpts.consensusSlack == 0.0) ? 1.0 : (1.0 - salmonOpts.consensusSlack);
  memCollector.setConsensusFraction(consensusFraction);
  memCollector.setHitFilterPolicy(salmonOpts.hitFilterPolicy);
  memCollector.setAltSkip(salmonOpts.mismatchSeedSkip);
  memCollector.setChainSubOptThresh(salmonOpts.pre_merge_chain_sub_thresh);

  //Initialize ksw aligner
  ksw2pp::KSW2Config config;
  config.dropoff = -1;
  config.gapo = salmonOpts.gapOpenPenalty;
  config.gape = salmonOpts.gapExtendPenalty;
  config.bandwidth = salmonOpts.dpBandwidth;
  config.flag = 0;
  config.flag |= KSW_EZ_RIGHT;
  config.flag |= KSW_EZ_SCORE_ONLY;
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
  aconf.allowOverhangSoftclip = salmonOpts.softclipOverhangs;
  aconf.allowSoftclip  = salmonOpts.softclip; 
  aconf.useAlignmentCache = !salmonOpts.disableAlignmentCache;
  aconf.alignmentMode = pufferfish::util::PuffAlignmentMode::SCORE_ONLY;

  // we actually care about the softclips in the cigar string 
  // if we are writing output and softclipping (or softclipping of overhangs) is enabled
  if ( (!salmonOpts.qmFileName.empty()) and (salmonOpts.softclip or salmonOpts.softclipOverhangs) ) {
    aconf.alignmentMode = pufferfish::util::PuffAlignmentMode::APPROXIMATE_CIGAR;
  }

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


inline void updateRefMappings(uint32_t tid, int32_t hitScore, bool isCompat, size_t idx,
                  const std::vector<Transcript>& transcripts,
                  int32_t invalidScore,
                  salmon::mapping_utils::MappingScoreInfo& msi) {
  auto& scores = msi.scores_;
  scores[idx] = hitScore;
  auto& t = transcripts[tid];
  bool isDecoy = t.isDecoy();
  double decoyCutoff = static_cast<int32_t>(msi.decoyThresh * msi.bestDecoyScore);

  //if (hitScore < decoyCutoff or (hitScore == invalidScore)) { }

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
  } else if ((hitScore > it->second.first) or (hitScore == it->second.first and isCompat)) {
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
                                       std::vector<pufferfish::util::JointMems>& jointHits,
                                       uint32_t readLen,
                                       uint32_t mateLen,
                                       bool singleEnd,
                                       bool tryAlign,
                                       bool hardFilter,
                                       double scoreExp,
                                       double minAlnProb,
                                       salmon::mapping_utils::MappingScoreInfo& msi,
                                       std::vector<pufferfish::util::QuasiAlignment>& jointAlignments) {

  auto invalidScore = std::numeric_limits<decltype(msi.bestDecoyScore)>::min();
  if (msi.bestDecoyScore == invalidScore) { msi.bestDecoyScore = invalidScore + 1 ;}

  int32_t decoyThreshold = static_cast<decltype(msi.bestDecoyScore)>(msi.decoyThresh * msi.bestDecoyScore);

  //auto filterScore = (bestDecoyScore < secondBestScore) ? secondBestScore : bestDecoyScore;
  auto& scores = msi.scores_;
  auto& perm = msi.perm_;

  // throw away any pairs for which we should not produce valid alignments :
  // ======
  // If we are doing soft-filtering (default), we remove those not exceeding the bestDecoyScore
  // If we are doing hard-filtering, we remove those less than the bestScore
  perm.erase(std::remove_if(perm.begin(), perm.end(),
                            [&scores, hardFilter, &msi, decoyThreshold, invalidScore](const std::pair<int32_t, int32_t>& idxtid) -> bool {
                              return !hardFilter ?  scores[idxtid.first] < decoyThreshold : scores[idxtid.first] < msi.bestScore;
                              //return !hardFilter ?  scores[idxtid.first] < filterScore : scores[idxtid.first] < bestScore; 
                            }), perm.end());

  // Unlike RapMap, pufferfish doesn't guarantee the hits computed above are in order
  // by transcript, so we find the permutation of indices that puts things in transcript
  // order.
  std::sort(perm.begin(), perm.end(), [](const std::pair<int32_t, int32_t>& p1,
                                         const std::pair<int32_t, int32_t>& p2) {
                                        return p1.second < p2.second;});

  // moving our alinged / score jointMEMs over to QuasiAlignment objects
  double bestScoreD = static_cast<double>(msi.bestScore);
  for (auto& idxTxp : perm) {
    int32_t ctr = idxTxp.first;
    int32_t tid = idxTxp.second;
    auto& jointHit = jointHits[ctr];

    double currScore = scores[ctr];
    double v = bestScoreD - currScore;
    // why -1?
    double estAlnProb = hardFilter ? -1.0 : std::exp(- scoreExp * v );
    // skip any alignment with aln prob < minAlnProb
    if (!hardFilter and (estAlnProb < minAlnProb)) { continue; }

    if (singleEnd or jointHit.isOrphan()) {
      readLen = jointHit.isLeftAvailable() ? readLen : mateLen;
      jointAlignments.emplace_back(tid,           // reference id
                                   jointHit.orphanClust()->getTrFirstHitPos(),     // reference pos
                                   jointHit.orphanClust()->isFw,     // fwd direction
                                   readLen, // read length
                                   jointHit.orphanClust()->cigar, // cigar string
                                   jointHit.fragmentLen,       // fragment length
                                   false);
      auto &qaln = jointAlignments.back();
      // NOTE : score should not be filled in from a double
      qaln.score = !tryAlign ? static_cast<int32_t >(jointHit.orphanClust()->coverage):jointHit.alignmentScore;
      qaln.estAlnProb(estAlnProb);
      // NOTE : wth is numHits?
      qaln.numHits = static_cast<uint32_t >(jointHits.size());//orphanClust()->coverage;
      qaln.mateStatus = jointHit.mateStatus;
      if(singleEnd) {
        qaln.mateLen = readLen;
        qaln.mateCigar.clear();
        qaln.matePos = 0;
        qaln.mateIsFwd = true;
        qaln.mateScore = 0;
        qaln.mateStatus = MateStatus::SINGLE_END;
      }
    } else {
      jointAlignments.emplace_back(tid,           // reference id
                                   jointHit.leftClust->getTrFirstHitPos(),     // reference pos
                                   jointHit.leftClust->isFw,     // fwd direction
                                   readLen, // read length
                                   jointHit.leftClust->cigar, // cigar string
                                   jointHit.fragmentLen,       // fragment length
                                   true);         // properly paired
      // Fill in the mate info
      auto &qaln = jointAlignments.back();
      qaln.mateLen = mateLen;
      qaln.mateCigar = jointHit.rightClust->cigar;
      qaln.matePos = static_cast<int32_t >(jointHit.rightClust->getTrFirstHitPos());
      qaln.mateIsFwd = jointHit.rightClust->isFw;
      qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
      // NOTE : wth is numHits?
      qaln.numHits = static_cast<uint32_t >(jointHits.size());
      // NOTE : score should not be filled in from a double
      qaln.score = !tryAlign ? static_cast<int32_t >(jointHit.leftClust->coverage):jointHit.alignmentScore;
      qaln.estAlnProb(estAlnProb);
      qaln.mateScore = !tryAlign ? static_cast<int32_t >(jointHit.rightClust->coverage):jointHit.mateAlignmentScore;
    }
  }
  // done moving our alinged / score jointMEMs over to QuasiAlignment objects
}


inline void filterAndCollectAlignmentsDecoy(
                                       std::vector<pufferfish::util::JointMems>& jointHits,
                                       uint32_t readLen,
                                       uint32_t mateLen,
                                       bool singleEnd,
                                       bool tryAlign,
                                       bool hardFilter,
                                       double scoreExp,
                                       double minAlnProb,
                                       salmon::mapping_utils::MappingScoreInfo& msi,
                                       std::vector<pufferfish::util::QuasiAlignment>& jointAlignments) {
// NOTE: this function should only be called in the case that we have valid decoy mappings to report.
// Currently, this happens only when there are *no valid non-decoy* mappings.
// Further, this function will only add equally *best* decoy mappings to the output jointAlignments object
// regardless of the the status of hardFilter (i.e. no sub-optimal decoy mappings will be reported).
(void) hardFilter;
(void) minAlnProb;
(void) scoreExp;
double estAlnProb = 1.0; //std::exp(-scoreExp * 0.0);
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
    qaln.score = !tryAlign ? static_cast<int32_t>(jointHit.leftClust->coverage)
                           : jointHit.alignmentScore;
    qaln.estAlnProb(estAlnProb);
    qaln.mateScore = !tryAlign
                         ? static_cast<int32_t>(jointHit.rightClust->coverage)
                         : jointHit.mateAlignmentScore;
  }
} // end for over best decoy hits
}


  } // namespace mapping_utils
} // namespace salmon

#endif //__SALMON_MAPPING_UTILS__
