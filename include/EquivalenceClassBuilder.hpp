#ifndef EQUIVALENCE_CLASS_BUILDER_HPP
#define EQUIVALENCE_CLASS_BUILDER_HPP

#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

// Logger includes
#include "spdlog/spdlog.h"

#include "SalmonUtils.hpp"
#include "TranscriptGroup.hpp"
#include "concurrentqueue.h"
#include "cuckoohash_map.hh"

struct TGValue {
  TGValue(const TGValue& o) {
    weights = o.weights;
    combinedWeights = o.combinedWeights;
    count = o.count;
    barcodeGroup = o.barcodeGroup;
  }

  TGValue(){}
  TGValue& operator=(const TGValue& o){
    weights = o.weights;
    combinedWeights = o.combinedWeights;
    count = o.count;
    //count.store(o.count.load());
    barcodeGroup = o.barcodeGroup;
    return *this;
  }

  TGValue(std::vector<double>& weightIn, uint64_t countIn)
      : weights(weightIn.begin(), weightIn.end()) {
    count = countIn;
  }

  //////////////////////////////////////////////////////////////////
  //constructor for handling barcodes
  TGValue(std::vector<double>& weightIn,
          uint64_t countIn, uint32_t barcode, uint64_t umi) :
    weights(weightIn.begin(), weightIn.end()) {
    count = countIn;
    barcodeGroup[barcode][umi] = 1;
  }
  //////////////////////////////////////////////////////////////////

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

  spp::sparse_hash_map<uint32_t,
                       spp::sparse_hash_map<uint64_t, uint32_t>> barcodeGroup;
};

class EquivalenceClassBuilder {
public:
  EquivalenceClassBuilder(std::shared_ptr<spdlog::logger> loggerIn)
      : logger_(loggerIn) {
    countMap_.reserve(1000000);
  }

  //~EquivalenceClassBuilder() {}

  void start() { active_ = true; }

  bool alv_finish(){
    active_ = false;
    size_t totalCount{0};
    auto lt = countMap_.lock_table();
    for (auto& kv : lt) {
      kv.second.normalizeAux();
      totalCount += kv.second.count;
    }

    logger_->info("Computed {} rich equivalence classes "
                  "for further processing", countMap_.size());
    logger_->info("Counted {} total reads in the equivalence classes ",
                  totalCount);
    return true;
  }

  bool finish() {
    active_ = false;
    size_t totalCount{0};
    auto lt = countMap_.lock_table();
    for (auto& kv : lt) {
      kv.second.normalizeAux();
      totalCount += kv.second.count;
      countVec_.push_back(kv);
    }

    logger_->info("Computed {} rich equivalence classes "
                  "for further processing",
                  countVec_.size());
    logger_->info("Counted {} total reads in the equivalence classes ",
                  totalCount);
    return true;
  }

  //////////////////////////////////////////////////////////////////
  //function for alevin barcode level count indexing
  inline void addBarcodeGroup(TranscriptGroup&& g,
                              std::vector<double>& weights,
                              uint32_t& barcode,
                              uint64_t& umi ){
    auto upfn = [&weights, &barcode, &umi](TGValue& x) -> void {
      // update the count
      x.count++;
      // update the weights
      for (size_t i = 0; i < x.weights.size(); ++i) {
        x.weights[i] += weights[i];
        //salmon::utils::incLoop(x.weights[i], weights[i]);
      }
      x.barcodeGroup[barcode][umi]++;
      /*
      try{
        x.barcodeGroup.at(barcode).at(umi)++;
      }
      catch (const std::out_of_range& e) {
        x.barcodeGroup[barcode][umi]  = 1;
        //cuckoo hash updates: deprecated
        //auto upbarfn = [](uint64_t &num) { ++num; };
        //barcodeGroup.upsert(barcode, upbarfn, 1);
      }
      */
    };

    // have to lock since tbb operator= is not concurrency safe
    TGValue v(weights, 1, barcode, umi);
    countMap_.upsert(g, upfn, v);
  }
  ////////////////////////////////////////////////////////////////

  inline void addGroup(TranscriptGroup&& g, std::vector<double>& weights) {

    auto upfn = [&weights](TGValue& x) -> void {
      // update the count
      x.count++;
      // update the weights
      for (size_t i = 0; i < x.weights.size(); ++i) {
        x.weights[i] += weights[i];
      }
    };
    TGValue v(weights, 1);
    countMap_.upsert(g, upfn, v);
  }

  cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher>& eqMap(){
    return countMap_;
  }

  std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec() {
    return countVec_;
  }

private:
  std::atomic<bool> active_;
  cuckoohash_map<TranscriptGroup, TGValue, TranscriptGroupHasher> countMap_;
  std::vector<std::pair<const TranscriptGroup, TGValue>> countVec_;
  std::shared_ptr<spdlog::logger> logger_;
};

#endif // EQUIVALENCE_CLASS_BUILDER_HPP

/** Unordered map implementation */
// std::unordered_map<TranscriptGroup, TGValue, TranscriptGroupHasher>
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
