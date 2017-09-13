#ifndef EFFECTIVE_LENGTH_STATS_HPP
#define EFFECTIVE_LENGTH_STATS_HPP

#include "Eigen/Dense"

/**
 * Record a weighted average of per-transcript observed effective lengths.
 * These can then be combined to compute an "expected" effective length.
 **/
class EffectiveLengthStats {
public:
  EffectiveLengthStats(size_t numTxps);
  void addFragment(uint32_t txID, uint32_t len, double logMass);
  double getExpectedEffectiveLength(uint32_t txID);
  Eigen::VectorXd getExpectedEffectiveLengths();
  void merge(const EffectiveLengthStats& other);
private:
  size_t numTxps_;
  //Eigen::Matrix<double, Eigen::Dynamic, 1> lengths_;
  Eigen::VectorXd lengths_;
  Eigen::VectorXd weights_;
};

#endif // EFFECTIVE_LENGTH_STATS_HPP
