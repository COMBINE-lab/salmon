/**
 * Model the bias in read start positions across a set of transcripts.
 * This class is similar to and inspired by the FragmentLengthDistribution
 * class, which was itself modified from lengthdistribution.cpp ---
 * originally written by Adam Roberts as part of the eXpress software.
 * Rob Patro; 2014, 2015
 */

#include "FragmentStartPositionDistribution.hpp"
#include "SalmonMath.hpp"
#include <boost/assign.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>

using namespace std;

FragmentStartPositionDistribution::FragmentStartPositionDistribution(
    uint32_t numBins)
    : numBins_(numBins), pmf_(numBins + 2), cmf_(numBins + 2),
      totMass_(salmon::math::LOG_1), isUpdated_(false), allowUpdates_(true),
      performingUpdate_(0) {

  using salmon::math::logAdd;
  double uniMass = log(1.0 / numBins_);
  pmf_[0] = salmon::math::LOG_0;
  cmf_[0] = salmon::math::LOG_0;
  for (size_t i = 1; i <= numBins_; ++i) {
    pmf_[i] = uniMass;
    cmf_[i] = log(static_cast<double>(i)) + uniMass;
  }
}

inline void logAddMass(std::atomic<double>& bin, double newMass) {
  double oldVal = bin.load();
  double newVal;
  do{
    newVal = salmon::math::logAdd(oldVal, newMass);
  } while (!bin.compare_exchange_strong(oldVal, newVal));
}

void FragmentStartPositionDistribution::addVal(int32_t hitPos, uint32_t txpLen,
                                               double mass) {

  ++performingUpdate_;
  {
    if (!allowUpdates_) {
      --performingUpdate_;
      return;
    }

    using salmon::math::logAdd;
    if (hitPos >= static_cast<int32_t>(txpLen)) {
      --performingUpdate_;
      return; // hit should happen within the transcript
    }

    if (hitPos < 0) {
      hitPos = 0;
    }
    double logLen = log(txpLen);

    // Modified from: https://github.com/deweylab/RSEM/blob/master/RSPD.h
    uint32_t i;
    // Fraction along the transcript where this hit occurs
    double a = hitPos * 1.0 / txpLen;
    double b;

    for (i = ((long long)hitPos) * numBins_ / txpLen + 1;
         i < (((long long)hitPos + 1) * numBins_ - 1) / txpLen + 1; i++) {
      b = i * 1.0 / numBins_;
      double updateMass = log(b - a) + logLen + mass;
      logAddMass(pmf_[i], updateMass);
      logAddMass(totMass_, updateMass);
      a = b;
    }
    b = (hitPos + 1.0) / txpLen;
    double updateMass = log(b - a) + logLen + mass;
    logAddMass(pmf_[i], updateMass);
    logAddMass(totMass_, updateMass);
  }
  --performingUpdate_;
}

double FragmentStartPositionDistribution::evalCDF(int32_t hitPos,
                                                  uint32_t txpLen) {
  int i = static_cast<int>((static_cast<double>(hitPos) * numBins_) / txpLen);
  double val = hitPos * (1.0 / txpLen) * numBins_;
  return (val - i < 1e-7)
             ? cmf_[i].load()
             : salmon::math::logAdd(cmf_[i], std::log(val - i) + pmf_[i + 1]);
}

void FragmentStartPositionDistribution::update() {
  if (isUpdated_) {
    return;
  }
  // TODO: Is this (thread)-safe yet?
  allowUpdates_ = false;
  std::lock_guard<std::mutex> lg(fspdMut_);
  // Make sure an update isn't being performed
  while (performingUpdate_) {
  }
  if (!isUpdated_) {
    for (uint32_t i = 1; i <= numBins_; i++) {
      pmf_[i] = pmf_[i] - totMass_;
      cmf_[i] = salmon::math::logAdd(cmf_[i - 1], pmf_[i]);
    }
    isUpdated_ = true;
  }
}

double FragmentStartPositionDistribution::
operator()(int32_t hitPos, uint32_t txpLen, double logEffLen) {
  if (hitPos < 0) {
    hitPos = 0;
  }
  uint32_t hitPosU = static_cast<uint32_t>(hitPos);
  assert(hitPosU < txpLen);
  if (hitPosU >= txpLen) {
    std::cerr << "\n\nhitPos = " << hitPosU << ", txpLen = " << txpLen
              << "!!\n\n\n";
    return salmon::math::LOG_0;
  }
  // If we haven't updated the CDF yet, then
  // just return log(1);
  if (!isUpdated_) {
    return -logEffLen;
  }

  double a = hitPosU * (1.0 / txpLen);

  double effLen = std::exp(logEffLen);
  if (effLen >= txpLen) {
    effLen = txpLen - 1;
  }

  double denom =
      evalCDF(static_cast<int32_t>(effLen), txpLen); // cmf_[numBins_];
  double cdfNext = evalCDF(hitPos + 1, txpLen);
  double cdfCurr = evalCDF(hitPos, txpLen);

  return ((denom >= salmon::math::LOG_EPSILON)
              ? salmon::math::logSub(cdfNext, cdfCurr) - denom
              : salmon::math::LOG_0);
}

bool FragmentStartPositionDistribution::logNumDenomMass(int32_t hitPos,
                                                        uint32_t txpLen,
                                                        double logEffLen,
                                                        double& logNum,
                                                        double& logDenom) {

  if (hitPos < 0) {
    hitPos = 0;
  }
  uint32_t uHitPos = static_cast<uint32_t>(hitPos);
  assert(uHitPos < txpLen);
  if (uHitPos >= txpLen) {
    std::cerr << "\n\nhitPos = " << hitPos << ", txpLen = " << txpLen
              << "!!\n\n\n";
    logNum = salmon::math::LOG_0;
    logDenom = salmon::math::LOG_0;
    return false;
  }
  // If we haven't updated the CDF yet, then
  // just return log(1);
  if (!isUpdated_) {
    logNum = std::log(1.0 / static_cast<double>(txpLen));
    logDenom = salmon::math::LOG_1;
    return true;
  }

  double effLen = std::exp(logEffLen);
  if (effLen >= txpLen) {
    effLen = txpLen - 1;
  }

  double denom = evalCDF(static_cast<int32_t>(effLen), txpLen);
  double cdfNext = evalCDF(hitPos + 1, txpLen);
  double cdfCurr = evalCDF(hitPos, txpLen);

  if (denom >= salmon::math::LOG_EPSILON) {
    logNum = salmon::math::logSub(cdfNext, cdfCurr);
    logDenom = denom;
    return true;
  } else {
    logNum = salmon::math::LOG_0;
    logDenom = salmon::math::LOG_0;
    return false;
  }
}

double FragmentStartPositionDistribution::totMass() const { return totMass_; }

std::string FragmentStartPositionDistribution::toString() const {
  std::stringstream ss;
  for (size_t i = 0; i < pmf_.size(); ++i) {
    ss << std::exp(pmf_[i] - totMass_);
    if (i != pmf_.size() - 1) {
      ss << '\t';
    }
  }
  ss << "\n";
  return ss.str();
}
