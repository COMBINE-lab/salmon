//
// Created by Fatemeh Almodaresi on 8/28/18.
//

#ifndef PUFFERFISH_CHAINFINDER_H
#define PUFFERFISH_CHAINFINDER_H

#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"
#include "CommonTypes.hpp"
#include "Util.hpp"

#include <sparsepp/spp.h>
#include "parallel_hashmap/phmap.h"
#include "metro/metrohash64.h"

class MemClusterer {
private:
  uint32_t maxAllowedRefsPerHit_ = 1000;
  double consensusFraction_ = 0.65;
  double pre_merge_chain_sub_thresh_{0.9};
  double inv_pre_merge_chain_sub_thresh_{1.0/0.9};

  using RefMemMap = pufferfish::util::CachedVectorMap<std::pair<pufferfish::common_types::ReferenceID, bool>, std::vector<pufferfish::util::MemInfo>, pufferfish::util::pair_hash>;

public:
  void set_chain_sub_opt_thresh_(double pre_merge_chain_sub_thresh, double inv_pre_merge_chain_sub_thresh);
  double chainSubOptThresh() const;
  void setConsensusFraction(double cf);
  double getConsensusFraction() const;
  void setMaxAllowedRefsPerHit(uint32_t max);
  uint32_t getMaxAllowedRefsPerHit();
  void setHitFilterPolicy(pufferfish::util::HitFilterPolicy hfp);
  pufferfish::util::HitFilterPolicy getHitFilterPolicy() const;

  size_t fillMemCollection(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                         //pufferfish::common_types::RefMemMapT& trMemMap,
                         RefMemMap& trMemMap,
                           std::vector<pufferfish::util::UniMemInfo> &memCollection, uint64_t firstDecoyIndex,
                           phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs);

  bool findOptChain(std::vector<std::pair<int, pufferfish::util::ProjectedHits>> &hits,
                    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>>& memClusters,
                    //phmap::flat_hash_map<pufferfish::common_types::ReferenceID, std::vector<pufferfish::util::MemCluster>> &memClusters,
                      uint32_t maxSpliceGap, std::vector<pufferfish::util::UniMemInfo> &memCollection, uint32_t readLen,
                    phmap::flat_hash_map<pufferfish::common_types::ReferenceID, bool>& other_end_refs,
                    bool hChain,
                    RefMemMap& trMemMap,
                    uint64_t firstDecoyIndex,
                    //pufferfish::common_types::RefMemMapT& trMemMap,
                    bool verbose = false);

private:
  itlib::small_vector<double> f;
  itlib::small_vector<int32_t> p;
  itlib::small_vector<uint8_t> keepMem;
  itlib::small_vector<uint64_t> memIndicesInReverse;
  itlib::small_vector<int32_t> bestChainEndList;
  itlib::small_vector<int32_t> chainQuerySig;
  pufferfish::util::HitFilterPolicy hitFilterPolicy_{pufferfish::util::HitFilterPolicy::FILTER_AFTER_CHAINING};
};

#endif //PUFFERFISH_CHAINFINDER_H
