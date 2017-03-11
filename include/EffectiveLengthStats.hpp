#ifndef EFFECTIVE_LENGTH_STATS_HPP
#define EFFECTIVE_LENGTH_STATS_HPP

#include "Eigen/Dense"

class EffectiveLengthStats {
public:
  EffectiveLengthStats(size_t numTxps);
  void addFragment(uint32_t txID, uint32_t len);
  double getExpectedEffectiveLength(uint32_t txID);
  Eigen::VectorXd getExpectedEffectiveLengths();
  void merge(const EffectiveLengthStats& other);
private:
  size_t numTxps_;
  Eigen::Matrix<uint64_t, Eigen::Dynamic, 1> lengths_;
  Eigen::VectorXi counts_;
};

#endif // EFFECTIVE_LENGTH_STATS_HPP
