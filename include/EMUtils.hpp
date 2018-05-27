#ifndef EM_UTILS_HPP
#define EM_UTILS_HPP

#include <vector>
#include "tbb/atomic.h"
#include "Transcript.hpp"

template <typename VecT>
void EMUpdate_(std::vector<std::vector<uint32_t>>& txpGroupLabels,
               std::vector<std::vector<double>>& txpGroupCombinedWeights,
               const std::vector<uint64_t>& txpGroupCounts,
               std::vector<Transcript>& transcripts, const VecT& alphaIn,
               VecT& alphaOut);

template <typename VecT>
double truncateCountVector(VecT& alphas, double cutoff);

#endif // EM_UTILS_HPP
