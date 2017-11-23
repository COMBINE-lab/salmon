#ifndef EFFECTIVE_LENGTH_STATS_HPP
#define EFFECTIVE_LENGTH_STATS_HPP

#include "Eigen/Dense"

/**
 * Record a weighted average of per-transcript observed effective lengths.
 * These can then be combined to compute an "expected" effective length.
 **/
class EffectiveLengthStats {
  using VectorXu = Eigen::Matrix<uint32_t, Eigen::Dynamic, 1>;

public:
  EffectiveLengthStats(size_t numTxps);
  void addFragment(uint32_t txID, uint32_t len, double logMass);
  uint32_t getObservedCount(uint32_t txID);
  double getExpectedEffectiveLength(uint32_t txID);
  Eigen::VectorXd getExpectedEffectiveLengths();
  void merge(const EffectiveLengthStats& other);

private:
  size_t numTxps_;
  Eigen::VectorXd lengths_;
  Eigen::VectorXd weights_;
  VectorXu counts_;
};

#endif // EFFECTIVE_LENGTH_STATS_HPP
