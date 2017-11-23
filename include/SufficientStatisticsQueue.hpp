#ifndef __SUFFICIENT_STATISTICS_QUEUE_HPP__
#define __SUFFICIENT_STATISTICS_QUEUE_HPP__

#include <cmath>
#include <mutex>
#include <queue>
#include <tuple>

#include "btree_map.h"

#include "SalmonMath.hpp"

class SufficientStatisticsQueue {

  using TranscriptID = uint32_t;

public:
  SufficientStatisticsQueue(uint32_t qlen = 10) : ct_(0), qlen_(qlen) {
    if (qlen_ > 0) {
      logL_ = std::log(static_cast<double>(qlen_));
    } else {
      std::cerr
          << "ERROR: cannot instantiate sufficient statistic queue of length 0";
    }
  }

  bool approxEq(double a, double b) { return std::abs(a - b) < 1e-3; }

  btree::btree_map<TranscriptID, std::tuple<double, size_t>>
  getSmoothedStats(btree::btree_map<TranscriptID, TranscriptAlignments>& Si) {
    using salmon::math::logAdd;
    using salmon::math::logSub;

    std::lock_guard<std::mutex> lg(ssmut_);

    ++ct_;
    if (ct_ % 200 == 0) {
      std::cerr << "\n\n sstatQueue_.size() = " << sstatQueue_.size()
                << ", suf stat map size = " << smoothedSufficientStats_.size()
                << "\n\n\n";
    }

    sstatQueue_.emplace();
    auto& currentSStat = sstatQueue_.back();

    // S_i^L = S_{i-1}^{L} + (S_{i} / L)
    for (auto kv = Si.begin(); kv != Si.end(); ++kv) {
      auto transcriptID = kv->first;
      auto& hits = kv->second;
      double contrib = hits.totalProb - logL_;
      currentSStat[transcriptID] = contrib;
      auto it = smoothedSufficientStats_.find(transcriptID);
      if (it == smoothedSufficientStats_.end()) {
        smoothedSufficientStats_[transcriptID] = std::make_tuple(contrib, 1);
      } else {
        std::get<0>(it->second) = logAdd(std::get<0>(it->second), contrib);
        std::get<1>(it->second) += 1;
      }
    }

    if (sstatQueue_.size() > qlen_) {
      // S_i^L -= (S_{i-L} / L)
      auto& Sil = sstatQueue_.front();
      for (auto kv = Sil.begin(); kv != Sil.end(); ++kv) {
        auto transcriptID = kv->first;
        auto& lprob = kv->second;
        auto it = smoothedSufficientStats_.find(transcriptID);
        if (it == smoothedSufficientStats_.end()) {
          // std::cerr << "\n\nWHAT\n\n\n";
        } else if (std::get<0>(it->second) <= lprob or
                   approxEq(
                       lprob,
                       std::get<0>(
                           it->second))) { // std::get<1>(it->second) == 1) {
          if (!approxEq(std::get<0>(it->second), lprob) and
              std::get<1>(it->second) > 1) {
            std::cerr << "hrmmmm .... remaining mass = "
                      << std::get<0>(it->second) << " < lprob = " << lprob
                      << "\n";
          }
          smoothedSufficientStats_.erase(it);
        } else {
          double before = std::get<0>(it->second);
          double update = logSub(std::get<0>(it->second), lprob);
          std::get<0>(it->second) = update;
          std::get<1>(it->second) -= 1;
          if (std::get<0>(it->second) == salmon::math::LOG_0) {
            smoothedSufficientStats_.erase(it);
          }
        }
      }
      sstatQueue_.pop();
    }

    // make sure this is a copy!
    return smoothedSufficientStats_;
    // lock is freed here
  }

private:
  size_t ct_;
  std::mutex ssmut_;
  std::queue<btree::btree_map<TranscriptID, double>> sstatQueue_;
  btree::btree_map<TranscriptID, std::tuple<double, size_t>>
      smoothedSufficientStats_;
  uint32_t qlen_;
  double logL_;
};

#endif // __SUFFICIENT_STATISTICS_QUEUE_HPP__
