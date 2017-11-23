#ifndef __FORGETTING_MASS_CALCULATOR__
#define __FORGETTING_MASS_CALCULATOR__

#include "SalmonMath.hpp"
#include "SalmonSpinLock.hpp"
#include "spdlog/spdlog.h"

class ForgettingMassCalculator {
public:
  ForgettingMassCalculator(double forgettingFactor = 0.65)
      : batchNum_(0), forgettingFactor_(forgettingFactor),
        logForgettingMass_(salmon::math::LOG_1), logForgettingMasses_({}),
        cumulativeLogForgettingMasses_({}) {}

  ForgettingMassCalculator(const ForgettingMassCalculator&) = delete;
  ForgettingMassCalculator(ForgettingMassCalculator&&) = default;
  ForgettingMassCalculator& operator=(ForgettingMassCalculator&&) = default;
  ForgettingMassCalculator& operator=(const ForgettingMassCalculator&) = delete;

  /** Precompute the log(forgetting mass) and cumulative log(forgetting mass)
   * for the first numMiniBatches batches / timesteps.
   */
  bool prefill(uint64_t numMiniBatches) {
    logForgettingMasses_.reserve(numMiniBatches);
    cumulativeLogForgettingMasses_.reserve(numMiniBatches);

    double fm = salmon::math::LOG_1;
    logForgettingMasses_.push_back(fm);
    cumulativeLogForgettingMasses_.push_back(fm);

    for (size_t i = 2; i < numMiniBatches - 1; ++i) {
      fm += forgettingFactor_ * std::log(static_cast<double>(i - 1)) -
            std::log(std::pow(static_cast<double>(i), forgettingFactor_) - 1);
      logForgettingMasses_.push_back(fm);
      // fill in cumulative mass
      cumulativeLogForgettingMasses_.push_back(
          salmon::math::logAdd(cumulativeLogForgettingMasses_.back(), fm));
    }
    return true;
  }

  double operator()() {
#if defined __APPLE__
    spin_lock::scoped_lock sl(ffMutex_);
#else
    std::lock_guard<std::mutex> lock(ffMutex_);
#endif
    ++batchNum_;
    if (batchNum_ > 1) {
      logForgettingMass_ +=
          forgettingFactor_ * std::log(static_cast<double>(batchNum_ - 1)) -
          std::log(std::pow(static_cast<double>(batchNum_), forgettingFactor_) -
                   1);
    }
    return logForgettingMass_;
  }

  /**
   *  Return the log(forgetting mass) and current timestep in the ouput
   *  variables logForgettingMass and currentMinibatchTimestep. If we
   *  haven't pre-computed the forgetting mass for the next timestep yet,
   *  then do it now.
   */
  void getLogMassAndTimestep(double& logForgettingMass,
                             uint64_t& currentMinibatchTimestep) {
#if defined __APPLE__
    spin_lock::scoped_lock sl(ffMutex_);
#else
    std::lock_guard<std::mutex> lock(ffMutex_);
#endif
    if (batchNum_ < logForgettingMasses_.size()) {
      currentMinibatchTimestep = batchNum_;
      logForgettingMass = logForgettingMasses_[batchNum_];
      ++batchNum_;
    } else {
      double fm =
          logForgettingMasses_.back() +
          forgettingFactor_ * std::log(static_cast<double>(batchNum_ - 1)) -
          std::log(std::pow(static_cast<double>(batchNum_), forgettingFactor_) -
                   1);

      logForgettingMasses_.push_back(fm);
      cumulativeLogForgettingMasses_.push_back(
          salmon::math::logAdd(cumulativeLogForgettingMasses_.back(), fm));

      currentMinibatchTimestep = batchNum_;
      logForgettingMass = logForgettingMasses_[batchNum_];
      ++batchNum_;
    }
  }

  // Retrieve the log(forgetting mass) at a particular timestep.  This
  // function assumes that the forgetting mass has already been computed
  // for this timestep --- otherwise, this will result in a fatal error.
  double logMassAt(uint64_t timestep) {
    if (timestep < logForgettingMasses_.size()) {
      return logForgettingMasses_[timestep];
    } else {
      spdlog::get("jointLog")
          ->error("Requested forgetting mass for timestep {} "
                  "where it has not yet been computed.  This "
                  "likely means that the ForgettingMassCalculator "
                  "class was being used incorrectly!  Please "
                  "report this crash on GitHub!\n",
                  timestep);
      std::exit(1);
      return salmon::math::LOG_0;
    }
  }

  // Retrieve the cumulative log(forgetting mass) at a particular timestep.
  // This function assumes that the forgetting mass has already been computed
  // for this timestep --- otherwise, this will result in a fatal error.
  double cumulativeLogMassAt(uint64_t timestep) {
    if (timestep < cumulativeLogForgettingMasses_.size()) {
      return cumulativeLogForgettingMasses_[timestep];
    } else {
      spdlog::get("jointLog")
          ->error("Requested cumulative forgetting mass for timestep {} "
                  "where it has not yet been computed.  This "
                  "likely means that the ForgettingMassCalculator "
                  "class was being used incorrectly!  Please "
                  "report this crash on GitHub!\n",
                  timestep);
      std::exit(1);
      return salmon::math::LOG_0;
    }
  }

  uint64_t getCurrentTimestep() { return batchNum_; }

private:
  uint64_t batchNum_;
  double forgettingFactor_;
  double logForgettingMass_;
  std::vector<double> logForgettingMasses_;
  std::vector<double> cumulativeLogForgettingMasses_;
#if defined __APPLE__
  spin_lock ffMutex_;
#else
  std::mutex ffMutex_;
#endif
};

#endif //__FORGETTING_MASS_CALCULATOR__
