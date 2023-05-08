#ifndef EQUIVALENCE_CLASS_BUILDER_HPP
#define EQUIVALENCE_CLASS_BUILDER_HPP

#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>
#include <atomic>

#include "cedar/TargetGroup.hpp"

struct TGValue {
  TGValue(const TGValue& o) {
    weights = o.weights;
    combinedWeights = o.combinedWeights;
    count = o.count;
  }

  TGValue(std::vector<double>& weightIn, uint64_t countIn)
      : weights(weightIn.begin(), weightIn.end()) {
    count = countIn;
  }

  // const is a lie
  void normalizeAux() const {
    double sumOfAux{0.0};
    for (size_t i = 0; i < weights.size(); ++i) {
      sumOfAux += weights[i];
    }
    double norm = 1.0 / sumOfAux;
    for (size_t i = 0; i < weights.size(); ++i) {
      weights[i] *= norm;
    }
  }

  mutable std::vector<double> weights;
  // The combined auxiliary and position weights.  These
  // are filled in by the inference algorithm.
  mutable std::vector<double> combinedWeights;
  uint64_t count{0};
};

class EquivalenceClassBuilder {
public:
  EquivalenceClassBuilder()//std::shared_ptr<spdlog::logger> loggerIn)
      //: logger_(loggerIn) {
  {
    countMap_.reserve(1000000);
  }

  //~EquivalenceClassBuilder() {}

  void start() { active_ = true; }

  bool finish() {
    active_ = false;
    size_t totalCount{0};
    countVec_.reserve(countMap_.size());
    //auto lt = countMap_.lock_table();
    for (auto& kv : countMap_) {
      kv.second.normalizeAux();
      totalCount += kv.second.count;
      countVec_.push_back(kv);
    }

    /*
    logger_->info("Computed {} rich equivalence classes "
                  "for further processing",
                  countVec_.size());
    logger_->info("Counted {} total reads in the equivalence classes ",
                  totalCount);
                  */
    return true;
  }

  inline void addGroup(TargetGroup&& g, std::vector<double>& weights) {
    auto it = countMap_.find(g);
    if (it != countMap_.end()) {
      it->second.count++;
      for (size_t i = 0; i < it->second.weights.size(); ++i) {
        it->second.weights[i] += weights[i];
      }
    } else {
      TGValue v(weights, 1);
      countMap_.insert(std::make_pair(g,v));
    }
  }

  inline void mergeUnfinishedEQB(EquivalenceClassBuilder& eqb) {
    for (auto &kv : eqb.countMap_) {
      auto it = countMap_.find(kv.first);
      if (it != countMap_.end()) {
        it->second.count+=kv.second.count;
        for (size_t i = 0; i < it->second.weights.size(); ++i) {
          it->second.weights[i] += kv.second.weights[i];
        }
      } else {
        countMap_.insert(std::make_pair(kv.first,kv.second));
      }
    }
  }

  std::vector<std::pair<const TargetGroup, TGValue>>& eqVec() {
    return countVec_;
  }
  std::unordered_map<TargetGroup, TGValue, TargetGroupHasher>& eqMap() {
    return countMap_;
  } 
private:
  std::atomic<bool> active_;
  std::unordered_map<TargetGroup, TGValue, TargetGroupHasher> countMap_;
  std::vector<std::pair<const TargetGroup, TGValue>> countVec_;
  //std::shared_ptr<spdlog::logger> logger_;
};

#endif // EQUIVALENCE_CLASS_BUILDER_HPP

/** Unordered map implementation */
// spp::sparse_hash_map<TranscriptGroup, TGValue, TranscriptGroupHasher>
// countMap_;  std::mutex mapMut_;
/*
bool finish() {
    // unordered_map implementation
    for (auto& kv : countMap_) {
        kv.second.normalizeAux();
        countVec_.push_back(kv);
    }
    return true;
}
*/

/*
inline void addGroup(TranscriptGroup&& g,
        std::vector<double>& weights) {

    // unordered_map implementation
    std::lock_guard<std::mutex> lock(mapMut_);
    auto it = countMap_.find(g);
    if (it == countMap_.end()) {
        TGValue v(weights, 1);
        countMap_.emplace(g, v);
    } else {
        auto& x = it->second;
        x.count++;
        for (size_t i = 0; i < x.weights.size(); ++i) {
            x.weights[i] =
                salmon::math::logAdd(x.weights[i], weights[i]);
        }
    }
}
*/
