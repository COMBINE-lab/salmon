#ifndef EM_UTILS_HPP
#define EM_UTILS_HPP

#include <vector>
#include <atomic>
#include "Transcript.hpp"

template <typename VecT>
void EMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               const std::vector<uint64_t>& txpGroupCounts,
               const VecT& alphaIn,
               VecT& alphaOut);

template <typename VecT>
void EMUpdate_augmented(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               const std::vector<uint64_t>& txpGroupCounts,
               const VecT& alphaIn,
               VecT& alphaOut,
               std::vector<uint32_t>& sampled_txps_counts,
               std::vector<bool>& eqClass_augmentable);
/**
 * set entries with values <= cutoff to 0.
 **/
template <typename VecT>
double truncateCountVector(VecT& alphas, double cutoff);

#endif // EM_UTILS_HPP
