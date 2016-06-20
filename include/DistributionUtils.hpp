#ifndef __DISTRIBUTION_UTILS__
#define __DISTRIBUTION_UTILS__

#include <vector>
#include <cstddef>
#include <cstdint>

// Do things involving distributions
class FragmentLengthDistribution;
class Transcript;

namespace distribution_utils {
  enum class DistributionSpace : uint8_t { LOG=0, LINEAR=1 };

  std::vector<double> correctionFactorsFromMass(std::vector<double>& mass, DistributionSpace inputSpace);
  void computeSmoothedEffectiveLengths(size_t maxLength,
				       std::vector<Transcript>& transcripts,
				       std::vector<double>& correctionFactors,
				       DistributionSpace outputSpace);
  std::vector<int32_t> samplesFromLogPMF(FragmentLengthDistribution* fld, int32_t numSamples);
}

#endif // __DISTRIBUTION_UTILS__
