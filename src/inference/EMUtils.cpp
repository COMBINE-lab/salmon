#include "EMUtils.hpp"
#include "SalmonUtils.hpp"

/**
 * Single-threaded EM-update routine for use in bootstrapping
 */
template <typename VecT>
void EMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               const std::vector<uint64_t>& txpGroupCounts,
               const VecT& alphaIn,
               VecT& alphaOut) {

  assert(alphaIn.size() == alphaOut.size());

  size_t numEqClasses = txpGroupLabels.size();
  for (size_t eqID = 0; eqID < numEqClasses; ++eqID) {
    uint64_t count = txpGroupCounts[eqID];
    // for each transcript in this class
    const std::vector<uint32_t>& txps = txpGroupLabels[eqID];
    const auto& auxs = txpGroupCombinedWeights[eqID];

    double denom = 0.0;
    size_t groupSize = txpGroupCombinedWeights[eqID].size(); // txps.size();
    // If this is a single-transcript group,
    // then it gets the full count.  Otherwise,
    // update according to our VBEM rule.
    if (BOOST_LIKELY(groupSize > 1)) {
      for (size_t i = 0; i < groupSize; ++i) {
        auto tid = txps[i];
        auto aux = auxs[i];
        double v = alphaIn[tid] * aux;
        denom += v;
      }

      if (denom <=  std::numeric_limits<double>::denorm_min()) {
        // tgroup.setValid(false);
      } else {
        double invDenom = count / denom;
        for (size_t i = 0; i < groupSize; ++i) {
          auto tid = txps[i];
          auto aux = auxs[i];
          double v = alphaIn[tid] * aux;
          if (!std::isnan(v)) {
            salmon::utils::incLoop(alphaOut[tid], v * invDenom);
          }
        }
      }
    } else {
      salmon::utils::incLoop(alphaOut[txps.front()], count);
    }
  }
}

template <typename VecT>
double truncateCountVector(VecT& alphas, double cutoff) {
  // Truncate tiny expression values
  double alphaSum = 0.0;

  for (size_t i = 0; i < alphas.size(); ++i) {
    if (alphas[i] <= cutoff) {
      alphas[i] = 0.0;
    }
    alphaSum += alphas[i];
  }
  return alphaSum;
}

template
void EMUpdate_<std::vector<double>>(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                                    std::vector<std::vector<double>>& txpGroupCombinedWeights,
                                    const std::vector<uint64_t>& txpGroupCounts,
                                    const std::vector<double>& alphaIn,
                                    std::vector<double>& alphaOut);

template
void EMUpdate_<std::vector<std::atomic<double>>>(std::vector<std::vector<uint32_t>>& txpGroupLabels,
                                                 std::vector<std::vector<double>>& txpGroupCombinedWeights,
                                                 const std::vector<uint64_t>& txpGroupCounts,
                                                 const std::vector<std::atomic<double>>& alphaIn,
                                                 std::vector<std::atomic<double>>& alphaOut);

template
double truncateCountVector<std::vector<double>>(std::vector<double>& alphas, double cutoff);

template
double truncateCountVector<std::vector<std::atomic<double>>>(std::vector<std::atomic<double>>& alphas, double cutoff);
