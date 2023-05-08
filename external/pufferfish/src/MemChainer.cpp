#include <numeric>
#include "MemChainer.hpp"
#include "itlib/small_vector.hpp"

// from "fastapprox"
// https://github.com/romeric/fastapprox/blob/master/fastapprox/src/fastlog.h
static inline float fastlog2(float x) {
  union {
      float f;
      uint32_t i;
  } vx = {x};
  union {
      uint32_t i;
      float f;
  } mx = {(vx.i & 0x007FFFFF) | 0x3f000000};
  float y = vx.i;
  y *= 1.1920928955078125e-7f;

  return y - 124.22551499f
         - 1.498030302f * mx.f
         - 1.72587999f / (0.3520887068f + mx.f);
}

// from "fastapprox"
// https://github.com/romeric/fastapprox/blob/master/fastapprox/src/fastlog.h
static inline float fasterlog2(float x) {
  union {
      float f;
      uint32_t i;
  } vx = {x};
  float y = vx.i;
  y *= 1.1920928955078125e-7f;
  return y - 126.94269504f;
}

void MemClusterer::set_chain_sub_opt_thresh_(const double pre_merge_chain_sub_thresh, const double inv_pre_merge_chain_sub_thresh){
  pre_merge_chain_sub_thresh_ = pre_merge_chain_sub_thresh;
  inv_pre_merge_chain_sub_thresh_ = inv_pre_merge_chain_sub_thresh;
}

double MemClusterer::chainSubOptThresh() const {
  return pre_merge_chain_sub_thresh_;
}

void MemClusterer::setConsensusFraction(double cf) { consensusFraction_ = cf; }
double MemClusterer::getConsensusFraction() const { return consensusFraction_; }

void MemClusterer::setHitFilterPolicy(pufferfish::util::HitFilterPolicy hfp) {
  hitFilterPolicy_ = hfp;
}

pufferfish::util::HitFilterPolicy MemClusterer::getHitFilterPolicy() const {
  return hitFilterPolicy_;
}

void MemClusterer::setMaxAllowedRefsPerHit(uint32_t maxh){
  maxAllowedRefsPerHit_ = maxh;
}

uint32_t MemClusterer::getMaxAllowedRefsPerHit() {
  return maxAllowedRefsPerHit_;
}

size_t MemClusterer::fillMemCollection(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                                     //pufferfish::common_types::RefMemMapT &trMemMap,
                                     RefMemMap& trMemMap,
                                     std::vector<pufferfish::util::UniMemInfo> &memCollection, uint64_t firstDecoyIndex,
                                     phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool> & /*other_end_refs*/) {
  using namespace pufferfish::common_types;
  if (hits.empty()) {
    return 0;
  }

  size_t maxNonDecoyHits{0};
  size_t totSize{0};
  for (auto &hit : core::range<decltype(hits.begin())>(hits.begin(), hits.end())) {
    auto &refs = hit.second.refRange;
    auto rs = refs.size();
    totSize +=  (static_cast<uint64_t>(rs) < maxAllowedRefsPerHit_) ? rs : 0;
  }

  // here we guarantee that even if later we fill up
  // every gap between the hits and before the first and after the last hit
  // still the memCollection size doesn't change and hence the pointers are valid
  //memCollection.reserve(maxAllowedRefsPerHit * 2 * hits.size() + 1);
  memCollection.reserve(totSize);

  for (auto &hit : core::range<decltype(hits.begin())>(hits.begin(), hits.end())) {
    auto &readPos = hit.first;
    auto &projHits = hit.second;
    // NOTE: here we rely on internal members of the ProjectedHit (i.e., member variables ending in "_").
    // Maybe we want to change the interface (make these members public or provide accessors)?
    auto &refs = projHits.refRange;
    if (static_cast<uint64_t>(refs.size()) < maxAllowedRefsPerHit_) {
      uint32_t mappings{0};
      memCollection.emplace_back(projHits.contigIdx_, projHits.contigOrientation_,
                                 readPos, projHits.k_, projHits.contigPos_,
                                 projHits.globalPos_ - projHits.contigPos_, projHits.contigLen_, pufferfish::util::ReadEnd::LEFT);
      auto memItr = std::prev(memCollection.end());
      for (auto &posIt : refs) {
      //If we want to let the the hits to the references also found by the other end to be accepted
      //if (static_cast<uint64_t>(refs.size()) < maxAllowedRefsPerHit or other_end_refs.find(posIt.transcript_id()) != other_end_refs.end() ) {
        const auto& refPosOri = projHits.decodeHit(posIt);
        auto tid = posIt.transcript_id();
        auto& refHits = trMemMap[std::make_pair(tid, refPosOri.isFW)];
        refHits.emplace_back(memItr, refPosOri.pos, refPosOri.isFW);
        auto nh = refHits.size();
        maxNonDecoyHits = (tid < firstDecoyIndex) ? std::max(nh, maxNonDecoyHits) : maxNonDecoyHits;
        mappings++;
      //}
      }
    }
  }
  return maxNonDecoyHits;
}

bool MemClusterer::findOptChain(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                                pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& memClusters,
                                //phmap::flat_hash_map<pufferfish::common_types::ReferenceID, std::vector<pufferfish::util::MemCluster>> &memClusters,
                                uint32_t maxSpliceGap, std::vector<pufferfish::util::UniMemInfo> &memCollection,
                                uint32_t readLen,
                                phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                                bool hChain,
                                RefMemMap& trMemMap,
                                uint64_t firstDecoyIndex,
                                //pufferfish::common_types::RefMemMapT& trMemMap,
                                bool /*verbose*/) {
  using namespace pufferfish::common_types;
  using pufferfish::util::HitFilterPolicy;
  //(void)verbose;
  //std::cerr << "pre_merge_chain_sub_thresh_ = " << pre_merge_chain_sub_thresh_ 
  //          << ", inv = " << inv_pre_merge_chain_sub_thresh_ << "\n\n";
  // Map from (reference id, orientation) pair to a cluster of MEMs.
  size_t maxHits = fillMemCollection(hits, trMemMap, memCollection, firstDecoyIndex, other_end_refs);
  if (maxHits == 0) {
    return false;
  }

  bool filterBefore = (hitFilterPolicy_ == HitFilterPolicy::FILTER_BEFORE_CHAINING) or
    (hitFilterPolicy_ == HitFilterPolicy::FILTER_BEFORE_AND_AFTER_CHAINING);
  bool filterAfter= (hitFilterPolicy_ == HitFilterPolicy::FILTER_AFTER_CHAINING) or
    (hitFilterPolicy_ == HitFilterPolicy::FILTER_BEFORE_AND_AFTER_CHAINING);


  double maxChainScore{0.0};
  int32_t signedReadLen = static_cast<int32_t>(readLen);
  for (auto hitIt = trMemMap.begin(); hitIt != trMemMap.end(); ++hitIt) {
    auto& trOri = hitIt->first;
    auto &tid = trOri.first;
    auto &isFw = trOri.second;
    auto &memList = *hitIt->second;
    size_t hits = memList.size();
    if (filterBefore and (hits < consensusFraction_ * maxHits)) { continue; }

    // sort memList according to mem reference positions
    std::sort(memList.begin(), memList.end(),
              [isFw](pufferfish::util::MemInfo &q1, pufferfish::util::MemInfo &q2) -> bool {
                  auto q1ref = q1.tpos + q1.extendedlen;
                  auto q2ref = q2.tpos + q2.extendedlen;
                  auto q1read = q1.rpos + q1.extendedlen;
                  auto q2read = q2.rpos + q2.extendedlen;
                  return q1ref != q2ref ? q1ref < q2ref :
                         (isFw ? q1read < q2read : q1read > q2read);// sort based on tpos
              });
    /*if (verbose) {

      std::cerr << "\ntid" << tid << " , isFw:" << isFw << "\n";
      for (auto &m : memList) {
        std::cerr << "\ttpos:" << m.tpos << " rpos:" << m.memInfo->rpos << " len:" << m.memInfo->memlen
                  << "\n";
      }
    }*/

    //auto minPosIt = memList.begin();
    // find the valid chains
    // Use variant of minimap2 scoring (Li 2018)
    // https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty191/4994778
    auto alpha = [](int32_t qdiff, int32_t rdiff, int32_t ilen) -> double {
        double score = ilen;
        double mindiff = (qdiff < rdiff) ? qdiff : rdiff;
        return (score < mindiff) ? score : mindiff;
    };

    auto beta = [maxSpliceGap](int32_t qdiff, int32_t rdiff, double avgseed) -> double {
        double l = qdiff - rdiff;
        uint32_t al = std::abs(l);
        if (qdiff <= 0 or rdiff <= 0 or (al > maxSpliceGap)) {
          return std::numeric_limits<double>::infinity();
        }
        // To penalize cases with organized gaps for reads such as
        // CTCCTCATCCTCCTCATCCTCCTCCTCCTCCTCCTCCTCCGCTGCCGCCGCCGACCGACTGAACCGCACCCGCCGCGCCGCACCGCCTCCAAGTCCCGGC
        // polyester simulated on human transcriptome. 0.01 -> 0.05
        return (l == 0) ? 0.0 : (0.05 * avgseed * al + 0.5 * fastlog2(static_cast<float>(al)));
    };

    constexpr const double bottomScore = std::numeric_limits<double>::lowest();
    double bestScore = bottomScore;
    int32_t bestChainEnd = -1;
    double avgseed = 31.0;
    f.clear();
    p.clear();
    keepMem.clear();
    bestChainEndList.clear();
    //auto lastHitId = static_cast<int32_t>(memList.size() - 1);

    // Compact mems before chaining.
    // UniMEMs can terminate because of the end of a contig, even if
    // there is still an exact match between the read and one or more references.
    // Here, we compact UniMEMs that would have constituted a larger (contiguous)
    // MEM with respect to the current reference.
    int32_t prev_qposi_end = 0;
    int32_t prev_rposi_end = 0;
    size_t currentMemIdx = 0;

    //int32_t prev_qposi_start = -1;
    //int32_t prev_rposi_start = -1;

    int32_t totLen{0};
    //int32_t lastStartPos{0};
    for (int32_t i = 0; i < static_cast<int32_t>(memList.size()); ++i) {
      auto &hi = memList[i];
      int32_t qposi_start = hi.isFw ? hi.rpos : signedReadLen - (hi.rpos + hi.extendedlen);
      int32_t rposi_start = hi.tpos;

      int32_t qposi_end = hi.isFw ? (hi.rpos + hi.extendedlen) : (signedReadLen - hi.rpos);
      int32_t rposi_end = hi.tpos + hi.extendedlen;

      int32_t overlap_read = (prev_qposi_end - qposi_start);
      int32_t overlap_ref = (prev_rposi_end - rposi_start);
      if (i > 0 and overlap_ref >= 0 and (overlap_ref == overlap_read)) {
        auto &lastMem = memList[currentMemIdx];
        uint32_t extension = rposi_end - prev_rposi_end;
        lastMem.extendedlen += extension;
        totLen += extension;
        if (!isFw) {
          lastMem.rpos = hi.rpos;
        }
        hi.extendedlen = std::numeric_limits<decltype(hi.extendedlen)>::max();
      } else {
        totLen += hi.extendedlen;
        currentMemIdx=i;
      }
      //prev_qposi_start = qposi_start;
      //prev_rposi_start = rposi_start;
      prev_qposi_end = qposi_end;
      prev_rposi_end = rposi_end;
    }

    //if (bestBinCoverage < maxBinCoverage * consensusFraction_) { continue; }
    //maxBinCoverage = std::max(bestBinCoverage, maxBinCoverage);

    memList.erase(std::remove_if(memList.begin(), memList.end(),
                                 [](pufferfish::util::MemInfo& m) {
                                   bool r = m.extendedlen == std::numeric_limits<decltype(m.extendedlen)>::max(); return r;
                                 }), memList.end());
    avgseed = totLen / static_cast<double>(memList.size());

    /*
    if (verbose) {
      std::cerr << "\ntid" << tid << " , isFw:" << isFw << "\n";
      for (auto &m : memList) {
        std::cerr << "\ttpos:" << m.tpos << " rpos:" << m.rpos << " len:" << m.extendedlen
                  << "\n";
      }
    }
    */

    p.reserve(memList.size());
    f.reserve(memList.size());
    for (int32_t i = 0; i < static_cast<int32_t>(memList.size()); ++i) {
      auto &hi = memList[i];

      int32_t qposi = hi.rpos + hi.extendedlen;
      int32_t rposi = hi.tpos + hi.extendedlen;

      double baseScore = static_cast<double>(hi.extendedlen);
      p.push_back(i);
      f.push_back(baseScore);

      // possible predecessors in the chain
      int32_t numRounds{2};
      (void) numRounds;
      for (int32_t j = i - 1; j >= 0; --j) {
        auto &hj = memList[j];

        int32_t qposj = hj.rpos + hj.extendedlen;
        int32_t rposj = hj.tpos + hj.extendedlen;

        int32_t qdiff = isFw ? qposi - qposj :
                        (qposj - hj.extendedlen) - (qposi - hi.extendedlen);
        int32_t rdiff = rposi - rposj;

        auto extensionScore = f[j] + alpha(qdiff, rdiff, hi.extendedlen) - beta(qdiff, rdiff, avgseed);

        bool extendWithJ = (extensionScore > f[i]);
        p[i] = extendWithJ ? j : p[i];
        f[i] = extendWithJ ? extensionScore : f[i];

        // HEURISTIC : if we connected this match to an earlier one
        // i.e. if we extended the chain.
        // This implements Heng Li's heuristic ---
        // "
        // We note that if anchor i is chained to j, chaining i to a predecessor of j
        // is likely to yield a lower score.
        // "
        // here we take this to the extreme, and stop at the first j to which we chain.
        // we can add a parameter "h" as in the minimap paper.  But here we expect the
        // chains of matches in short reads to be short enough that this may not be worth it.
        if (hChain and p[i] < i) {
          numRounds--;
          if (numRounds <= 0) { break; }
        }
        // If the last two hits are too far from each other, we are sure that 
        // every other hit will be even further since the mems are sorted
        if (rdiff > signedReadLen * 2) {
          break;
        }
        // Mohsen: This heuristic hurts the accuracy of the chain in the case of this read:
        // TGAACGCTCTATGATGTCAGCCTACGAGCGCTCTATGATGTTAGCCTACGAGCGCTCTATGATGTCCCCTATGGCTGAGCGCTCTATGATGTCAGCTTAT
        // from Polyester simalted sample aligning to the human transcriptome
      }
      
      if (f[i] > bestScore * inv_pre_merge_chain_sub_thresh_) {
        bestScore = f[i];
        bestChainEnd = i;
        bestChainEndList.clear();
        bestChainEndList.push_back(bestChainEnd);
      } else if (f[i] >= bestScore * pre_merge_chain_sub_thresh_) {
        bestChainEndList.push_back(i);
        bestScore = std::max(f[i], bestScore);
      }
    }

    // because we were always comparing to the _current_ best chaining score, go through
    // and remove chains that may have passed at the time, but should not have passed globally
    bestChainEndList.erase(std::remove_if(bestChainEndList.begin(), bestChainEndList.end(),
      [this, bestScore](const int32_t& chain_ind) -> bool { 
        return this->f[chain_ind] < this->pre_merge_chain_sub_thresh_ * bestScore; 
      }), bestChainEndList.end());
    

    // early exit if this doesn't seem a promising chain
    if (filterAfter and bestScore < maxChainScore * consensusFraction_) { continue; }
    maxChainScore = std::max(bestScore, maxChainScore);

    // Do backtracking
    itlib::small_vector<uint8_t> seen(f.size(), 0);
    for (auto bestChainEnd : bestChainEndList) {
      if (bestChainEnd >= 0) {
        bool shouldBeAdded = true;
        memIndicesInReverse.clear();
        auto lastPtr = p[bestChainEnd];
        
        // keep track of the actual end index here, because 
        // we will need it to record the correct score later,
        // but we modify it when walking the chain.
        const auto currentChainEndInd = bestChainEnd;

        while (lastPtr < bestChainEnd) {
          if (seen[bestChainEnd] > 0) {
            shouldBeAdded = false;
            //break;
          }
          memIndicesInReverse.push_back(bestChainEnd);
          seen[bestChainEnd] = 1;
          bestChainEnd = lastPtr;
          lastPtr = p[bestChainEnd];
          //lastPtr = bestChainEnd;
          //bestChainEnd = p[bestChainEnd];
        }
        if (seen[bestChainEnd] > 0) {
          shouldBeAdded = false;
        }
        memIndicesInReverse.push_back(bestChainEnd);
        if (shouldBeAdded) {
          // @fataltes --- is there a reason we were inserting here before rather than
          // pushing back?
          memClusters[tid].push_back(pufferfish::util::MemCluster(isFw, signedReadLen));
          auto& justAddedCluster = memClusters[tid].back();
          for (auto it = memIndicesInReverse.rbegin(); it != memIndicesInReverse.rend(); it++) {
            justAddedCluster.addMem(memList[*it].memInfo, memList[*it].tpos,
                                       memList[*it].extendedlen, memList[*it].rpos, isFw);
          }
          justAddedCluster.coverage = f[currentChainEndInd];
          if (justAddedCluster.coverage == signedReadLen) {
            justAddedCluster.perfectChain = true;
          }
          /*
          if (verbose)
            std::cerr<<"Added position: " << memClusters[tid][0].coverage << " " << memClusters[tid][0].mems[0].tpos << "\n";
          */
        }
        //minPosIt += lastPtr;
      } else {
        // should not happen
        std::cerr << "[FATAL] : Cannot find any valid chain for quasi-mapping\n";
        std::cerr << "num hits = " << memList.size() << "\n";
        std::cerr << "bestChainEnd = " << bestChainEnd << "\n";
        std::cerr << "bestChainScore = " << bestScore << "\n";
        std::exit(1);
      }
    }

    for (auto & memClust : memClusters[tid]) {
      auto &memList = memClust.mems;
      size_t nmem = 0; // the number of mems that will be in this chain
      for (int32_t i = 0; i < static_cast<int32_t>(memList.size()); ++i) {
        auto &hi = memList[i];
        //chainOfInterest = /*chainOfInterest or */(hi.rpos == 1 and hi.tpos == 163 and tid == 151214);
        int32_t qposi_start = hi.isFw ? hi.rpos : signedReadLen - (hi.rpos + hi.extendedlen);
        int32_t rposi_start = hi.tpos;

        int32_t qposi_end = hi.isFw ? (hi.rpos + hi.extendedlen) : (signedReadLen - hi.rpos);
        int32_t rposi_end = hi.tpos + hi.extendedlen;

        int32_t overlap_read = (prev_qposi_end - qposi_start);
        int32_t overlap_ref = (prev_rposi_end - rposi_start);
        if (i > 0 and overlap_ref >= 0 and (overlap_ref == overlap_read)) {
          auto &lastMem = memList[currentMemIdx];
          uint32_t extension = rposi_end - prev_rposi_end;
          lastMem.extendedlen += extension;
          if (!isFw) {
            lastMem.rpos = hi.rpos;
          }
          hi.extendedlen = std::numeric_limits<decltype(hi.extendedlen)>::max();
        } else {
          ++nmem;
          currentMemIdx=i;
        }
        //prev_qposi_start = qposi_start;
        //prev_rposi_start = rposi_start;
        prev_qposi_end = qposi_end;
        prev_rposi_end = rposi_end;
      }

      chainQuerySig.clear();
      chainQuerySig.reserve(2*nmem);
      int32_t prev_tpos = -1;
      memList.erase(std::remove_if(memList.begin(), memList.end(),
                                   [this, signedReadLen, &prev_tpos](pufferfish::util::MemInfo& m) {
                                       bool r = m.extendedlen == std::numeric_limits<decltype(m.extendedlen)>::max(); 
                                       if (!r) {
                                        int32_t qposi_start = m.isFw ? m.rpos : signedReadLen - (m.rpos + m.extendedlen);
                                        int32_t rdiff = (prev_tpos == -1) ? 0 : (m.tpos - prev_tpos);
                                        prev_tpos = m.tpos;
                                        this->chainQuerySig.push_back(qposi_start);
                                        this->chainQuerySig.push_back(rdiff);
                                       }
                                       return r;
                                   }), memList.end());

      MetroHash64::Hash(reinterpret_cast<uint8_t *>(const_cast<int32_t*>(chainQuerySig.data())), 
                        chainQuerySig.size() * sizeof(int32_t), 
                        reinterpret_cast<uint8_t *>(&memClust.queryChainHash), 0);
    }

  }

  // since the maxChainScore was changing while we were computing the chains
  // do one final pass over the chains we collected to remove any that
  // should not be considered according to the final threshold.
  if (filterAfter) {
    auto cfrac = consensusFraction_;
    for (auto& mc : memClusters) {
      auto* clusters_for_txp = mc.second;
      clusters_for_txp->erase(
          std::remove_if(clusters_for_txp->begin(), clusters_for_txp->end(),
                      [maxChainScore, cfrac](const pufferfish::util::MemCluster& c) {
                        return c.coverage < maxChainScore * cfrac;
                      }),
          clusters_for_txp->end());
    }
  }

  return true;
}
