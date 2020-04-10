#ifndef EQUIVALENCE_CLASS_BUILDER_HPP
#define EQUIVALENCE_CLASS_BUILDER_HPP

#include <memory>
#include <mutex>
#include <thread>
#include <unordered_map>
#include <vector>

// Logger includes
#include "spdlog/spdlog.h"
#include "nonstd/optional.hpp"

#include "SalmonUtils.hpp"
#include "TranscriptGroup.hpp"
#include "concurrentqueue.h"
#include "cuckoohash_map.hh"
#include "pufferfish/sparsepp/spp.h"

struct EmptyBarcodeMapType {};
using SparseBarcodeMapType = spp::sparse_hash_map<uint32_t, spp::sparse_hash_map<uint64_t, uint32_t>>;
using BarcodeT = uint32_t;
using UMIT = uint64_t;

/**
 * NOTE : think of a potentially safer implementation of the barcode / non-barcode
 * version here, like using CRTP.
 **/
struct SCTGValue {
  SCTGValue(const SCTGValue& o) {
    count = o.count;
    barcodeGroup = o.barcodeGroup;
  }

  SCTGValue(){}
  SCTGValue& operator=(const SCTGValue& o){
    count = o.count;
    barcodeGroup = o.barcodeGroup;
    return *this;
  }

  SCTGValue(uint64_t countIn)
    : count(countIn) {}

  SCTGValue(std::vector<double>&, int)
  { std::cerr<<"invalid initialization in eqbuilder"<<std::endl; exit(1); }
  //////////////////////////////////////////////////////////////////
  //constructor for handling barcodes
  SCTGValue(uint64_t countIn, uint32_t barcode, uint64_t umi) {
    count = countIn;
    barcodeGroup[barcode][umi] = 1;
  }

  SCTGValue(uint32_t countIn, uint32_t barcode,
            uint64_t umi, bool /*bulkUpdate*/) {
    count = countIn;
    barcodeGroup[barcode][umi] = countIn;
  }

  void updateBarcodeGroup(BarcodeT barcode, UMIT umi) {
    barcodeGroup[barcode][umi]++;
  }

  void updateBarcodeGroup(BarcodeT barcode, UMIT umi, uint32_t count) {
    barcodeGroup[barcode][umi] += count;
  }
  //////////////////////////////////////////////////////////////////

  // const is a lie
  void normalizeAux() const {}

  uint64_t count{0};
  SparseBarcodeMapType barcodeGroup;
};


struct TGValue {
  TGValue(const TGValue& o) {
    weights = o.weights;
    combinedWeights = o.combinedWeights;
    count = o.count;
  }

  TGValue(){}
  TGValue& operator=(const TGValue& o){
    weights = o.weights;
    combinedWeights = o.combinedWeights;
    count = o.count;
    //count.store(o.count.load());
    return *this;
  }

  TGValue(std::vector<double>& weightIn, uint64_t countIn)
      : weights(weightIn.begin(), weightIn.end()) {
    count = countIn;
  }

  //////////////////////////////////////////////////////////////////
  //constructor for handling barcodes
  TGValue(std::vector<double>& weightIn,
          uint64_t countIn, uint32_t /*barcode*/, uint64_t /*umi*/) :
    weights(weightIn.begin(), weightIn.end()) {
    count = countIn;
  }
  //////////////////////////////////////////////////////////////////

  // We need this because otherwise the template will complain ... this **could be**
  // be instantiated, but isn't.  Figure out a cleaner way to do this;
  void updateBarcodeGroup(BarcodeT /*bc*/, UMIT /*umi*/) {}
  TGValue(int, BarcodeT /*bc*/, UMIT /*umi*/)
  { std::cerr<<"invalid initialization in eqbuilder"<<std::endl; exit(1); }

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

template <typename TGValueType>
class EquivalenceClassBuilder {
public:
  EquivalenceClassBuilder(std::shared_ptr<spdlog::logger> loggerIn, uint32_t maxResizeThreads)
      : logger_(loggerIn) {
    countMap_.max_num_worker_threads(maxResizeThreads);
    countMap_.reserve(1000000);
  }

  //~EquivalenceClassBuilder() {}
  void setMaxResizeThreads(uint32_t t) { countMap_.max_num_worker_threads(t); }
  uint32_t getMaxResizeThreads() const { return countMap_.max_num_worker_threads(); }

  void start() { active_ = true; }

  bool alv_finish(){
    active_ = false;
    size_t totalCount{0};
    auto lt = countMap_.lock_table();
    for (auto& kv : lt) {
      totalCount += kv.second.count;
    }

    logger_->info("Computed {:n} rich equivalence classes "
                  "for further processing", countMap_.size());
    logger_->info("Counted {:n} total reads in the equivalence classes ",
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

    logger_->info("Computed {:n} rich equivalence classes "
                  "for further processing",
                  countVec_.size());
    logger_->info("Counted {:n} total reads in the equivalence classes ",
                  totalCount);
    return true;
  }

  //////////////////////////////////////////////////////////////////
  //function for alevin barcode level count indexing
  inline void addBarcodeGroup(TranscriptGroup&& g,
                              uint32_t& barcode,
                              uint64_t& umi ){
    auto upfn = [&barcode, &umi](TGValueType& x) -> void {
      // update the count
      x.count++;
      // update the weights
      x.updateBarcodeGroup(barcode, umi);
    };

    // have to lock since tbb operator= is not concurrency safe
    TGValueType v(1, barcode, umi);
    countMap_.upsert(g, upfn, v);
  }
  ////////////////////////////////////////////////////////////////


  // If we have range-factorized equivalence classes, then
  // TranscriptGroup.txps.size() is *not* equal to the number of transcripts in
  // this equivalence class.  This function provides a generic way to get the
  // actual number of transcripts that label each equivalence class.
  // NOTE:  It is only valid to call this function once the finish() method has
  // been called on the EquivalenceClassBuilder.
  inline size_t getNumTranscriptsForClass(size_t eqIdx) const;

  inline void addGroup(TranscriptGroup&& g, std::vector<double>& weights);

  inline void populateTargets(std::vector<std::vector<uint32_t>>& eqclasses,
                              std::vector<std::vector<double>>& auxs_vals,
                              std::vector<uint32_t>& eqclass_counts,
                              std::vector<Transcript>& transcripts);

  libcuckoo::cuckoohash_map<TranscriptGroup, TGValueType, TranscriptGroupHasher>& eqMap(){
    return countMap_;
  }

  std::vector<std::pair<const TranscriptGroup, TGValueType>>& eqVec() {
    return countVec_;
  }

  // The returned value is only valid when the class is finalized!
  nonstd::optional<size_t> numEqClasses() const {
    return (active_) ? nonstd::nullopt : nonstd::optional<size_t>(countVec_.size());
  }

private:
  std::atomic<bool> active_;
  libcuckoo::cuckoohash_map<TranscriptGroup, TGValueType, TranscriptGroupHasher> countMap_;
  std::vector<std::pair<const TranscriptGroup, TGValueType>> countVec_;
  std::shared_ptr<spdlog::logger> logger_;
};

template <>
inline void EquivalenceClassBuilder<TGValue>::addGroup(TranscriptGroup&& g,
                                                       std::vector<double>& weights) {
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

template <>
inline void EquivalenceClassBuilder<TGValue>::populateTargets(
                                      std::vector<std::vector<uint32_t>>& eqclasses,
                                      std::vector<std::vector<double>>& auxs_vals,
                                      std::vector<uint32_t>& eqclass_counts,
                                      std::vector<Transcript>& transcripts) {
  for (size_t i = 0; i < eqclass_counts.size(); ++i) {
    uint32_t count = eqclass_counts[i];
    auto& tids = eqclasses[i];

    TGValue val(auxs_vals[i], count);
    TranscriptGroup tgroup(tids);

    countVec_.emplace_back(std::make_pair(std::move(tgroup), val));
    for (uint32_t tid: tids) {
      transcripts[tid].addTotalCount(count);
    }

    if ( tids.size() == 1 ) { transcripts[tids[0]].addUniqueCount(count); }
  }
}

template <>
inline size_t EquivalenceClassBuilder<TGValue>::getNumTranscriptsForClass(size_t eqIdx) const {
  return countVec_[eqIdx].second.weights.size();
}

template <>
inline size_t EquivalenceClassBuilder<SCTGValue>::getNumTranscriptsForClass(size_t eqIdx) const {
  return countVec_[eqIdx].first.txps.size();
}

// explicit instantiations
template class EquivalenceClassBuilder<TGValue>;
template class EquivalenceClassBuilder<SCTGValue>;

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
