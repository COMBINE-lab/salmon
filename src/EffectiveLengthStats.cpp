#include "EffectiveLengthStats.hpp"

EffectiveLengthStats::EffectiveLengthStats(size_t numTxps) : numTxps_(numTxps), lengths_(numTxps), counts_(numTxps) {
  lengths_.setZero();
  counts_.setZero();
}

void EffectiveLengthStats::addFragment(uint32_t txID, uint32_t len) {
  lengths_(txID) += len;
  counts_(txID) += 1;
}

double EffectiveLengthStats::getExpectedEffectiveLength(uint32_t txID) {
  return (counts_(txID) > 0) ? static_cast<double>(lengths_(txID)) / counts_(txID) : 0.01;
}

Eigen::VectorXd EffectiveLengthStats::getExpectedEffectiveLengths() {
  Eigen::VectorXd eel(numTxps_);
  eel.setZero();
  for (size_t i = 0; i < numTxps_; ++i) { eel(i) = getExpectedEffectiveLength(i); }
  return eel;
}

void EffectiveLengthStats::merge(const EffectiveLengthStats& other) {
  lengths_ += other.lengths_;
  counts_ += other.counts_;
}
