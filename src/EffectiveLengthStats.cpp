#include "EffectiveLengthStats.hpp"

EffectiveLengthStats::EffectiveLengthStats(size_t numTxps) : numTxps_(numTxps), lengths_(numTxps), weights_(numTxps) {
  lengths_ = Eigen::VectorXd::Constant(salmon::math::LOG_0);
  weights_ = Eigen::VectorXd::Constant(salmon::math::LOG_0);
}

void EffectiveLengthStats::addFragment(uint32_t txID, uint32_t len, double logMass) {
  len = (len >= 1) ? len : 1;
  const double logLen = std::log(static_cast<double>(len));
  lengths_(txID) = salmon::math::logAdd(lengths_(txID), logLen);
  weights_(txID) = salmon::math::logAdd(weights_(txID), logMass);
}

double EffectiveLengthStats::getExpectedEffectiveLength(uint32_t txID) {
  return (!salmon::math::isLog0(weights_(txID))) ? std::exp(lengths_(txID) - weights_(txID)) : 0.01;
  //return (weights_(txID) > 0) ? static_cast<double>(lengths_(txID)) / weights_(txID) : 0.01;
}

Eigen::VectorXd EffectiveLengthStats::getExpectedEffectiveLengths() {
  // expected effective lengths
  Eigen::VectorXd eel(numTxps_);
  eel.setZero();
  for (size_t i = 0; i < numTxps_; ++i) { eel(i) = getExpectedEffectiveLength(i); }
  return eel;
}

void EffectiveLengthStats::merge(const EffectiveLengthStats& other) {
  for (size_t i = 0; i < weights_.size(); ++i) {
    lengths_(i) = salmon::math::logAdd(lengths_(i), other.lengths_(i));
    weights_(i) = salmon::math::logAdd(weights_(i), other.weights_(i));
  }
}
