#ifndef SALMON_INTERNAL_QUANT_BIAS_LIBRARY_STATE_HPP
#define SALMON_INTERNAL_QUANT_BIAS_LIBRARY_STATE_HPP

#include <array>
#include <atomic>
#include <memory>
#include <utility>
#include <vector>

#include "salmon/internal/model/FragmentStartPositionDistribution.hpp"
#include "salmon/internal/model/GCFragModel.hpp"
#include "salmon/internal/model/LibraryTypeDetector.hpp"
#include "salmon/internal/model/ReadKmerDist.hpp"
#include "salmon/internal/model/SBModel.hpp"
#include "salmon/internal/model/SimplePosBias.hpp"
#include "salmon/internal/util/DistributionUtils.hpp"
#include "salmon/internal/util/UtilityFunctions.hpp"

namespace salmon::quant {

template <typename LibraryStateT>
struct BiasLibraryState {
  BiasLibraryState(LibraryStateT libraryIn, size_t numCondBins,
                   size_t numFragGCBins, size_t numBiasBins = 5)
      : library(std::move(libraryIn)),
        fragmentStartDists(numBiasBins),
        posBiasFW(numBiasBins),
        posBiasRC(numBiasBins),
        posBiasExpectFW(numBiasBins),
        posBiasExpectRC(numBiasBins),
        observedGC(numCondBins, numFragGCBins,
                   distribution_utils::DistributionSpace::LOG),
        expectedGC(numCondBins, numFragGCBins,
                   distribution_utils::DistributionSpace::LOG),
        expectedSeqBias(constExprPow(4, readBias[0].getK()), 1.0) {}

  LibraryStateT library;
  std::vector<double> conditionalMeans;

  std::vector<FragmentStartPositionDistribution> fragmentStartDists;

  std::vector<SimplePosBias> posBiasFW;
  std::vector<SimplePosBias> posBiasRC;
  std::vector<SimplePosBias> posBiasExpectFW;
  std::vector<SimplePosBias> posBiasExpectRC;

  double gcFracFwd{-1.0};
  GCFragModel observedGC;
  GCFragModel expectedGC;

  std::array<ReadKmerDist<6, std::atomic<uint32_t>>, 2> readBias;
  std::array<SBModel, 2> readBiasModelObserved;
  std::array<SBModel, 2> readBiasModelExpected;
  std::vector<double> expectedSeqBias;

  std::unique_ptr<LibraryTypeDetector> detector{nullptr};
};

} // namespace salmon::quant

#endif
