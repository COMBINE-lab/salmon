#include "DistributionUtils.hpp"
#include "FragmentLengthDistribution.hpp"
#include "Transcript.hpp"

#include <random>

namespace distribution_utils {

std::vector<double> correctionFactorsFromMass(std::vector<double>& mass,
                                              DistributionSpace inputSpace) {
  auto maxLen = mass.size();

  std::vector<double> correctionFactors(maxLen, 0.0);
  std::vector<double> vals(maxLen, 0.0);
  std::vector<double> multiplicities(maxLen, 0);

  multiplicities[0] = mass[0];

  double v{0.0};
  for (size_t i = 1; i < maxLen; ++i) {
    v = mass[i];
    vals[i] = static_cast<double>(v * i) + vals[i - 1];
    multiplicities[i] = v + multiplicities[i - 1];
    if (multiplicities[i] > 0) {
      correctionFactors[i] = vals[i] / multiplicities[i];
    }
  }
  return correctionFactors;
}

void computeSmoothedEffectiveLengths(size_t maxLength,
                                     std::vector<Transcript>& transcripts,
                                     std::vector<double>& correctionFactors,
                                     DistributionSpace outputSpace) {

  auto maxLen = maxLength;

  for (auto& txp : transcripts) {
    auto origLen = static_cast<double>(txp.RefLength);
    double correctionFactor = (origLen >= maxLen)
                                  ? correctionFactors[maxLen - 1]
                                  : correctionFactors[origLen];

    double effLen = origLen - correctionFactor + 1.0;
    if (effLen < 1.0) {
      effLen = origLen;
    }

    if (outputSpace == DistributionSpace::LOG) {
      txp.setCachedLogEffectiveLength(std::log(effLen));
    } else {
      txp.EffectiveLength = effLen;
    }
  }
}

std::vector<int32_t> samplesFromLogPMF(FragmentLengthDistribution* fld,
                                       int32_t numSamples) {
  std::vector<double> logPMF;
  size_t minVal;
  size_t maxVal;
  double logFLDMean = fld->mean();
  fld->dumpPMF(logPMF, minVal, maxVal);
  double sum = salmon::math::LOG_0;
  for (auto v : logPMF) {
    sum = salmon::math::logAdd(sum, v);
  }
  for (auto& v : logPMF) {
    v -= sum;
  }

  // Create the non-logged pmf
  std::vector<double> pmf(maxVal + 1, 0.0);
  for (size_t i = minVal; i < maxVal; ++i) {
    pmf[i] = std::exp(logPMF[i - minVal]);
  }

  // generate samples
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<int32_t> dist(pmf.begin(), pmf.end());

  std::vector<int32_t> samples(pmf.size());
  for (int32_t i = 0; i < numSamples; ++i) {
    ++samples[dist(gen)];
  }
  return samples;
}
} // namespace distribution_utils
