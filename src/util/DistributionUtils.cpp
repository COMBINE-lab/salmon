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

    double effLen = origLen - correctionFactor;
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

DistSummary samplesFromLogPMF(FragmentLengthDistribution* fld,
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
  double mean = salmon::math::LOG_0;
  double var = 0.0;
  std::vector<double> pmf(maxVal + 1, 0.0);
  for (size_t i = minVal; i < maxVal; ++i) {
    mean = (i > 0) ? (salmon::math::logAdd(mean, std::log(i)+logPMF[i-minVal])) : mean;
    pmf[i] = std::exp(logPMF[i - minVal]);
    var += pmf[i] * (i*i);
  }
  mean = std::exp(mean);
  var -= mean*mean;
  double sd = std::sqrt(var);

  // generate samples
  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    std::random_device rd("/dev/urandom");
  #else
    std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 2020012

  std::mt19937 gen(rd());
  std::discrete_distribution<int32_t> dist(pmf.begin(), pmf.end());

  DistSummary ds(mean, sd, pmf.size());
  std::vector<int32_t> samples(pmf.size());
  for (int32_t i = 0; i < numSamples; ++i) {
    ++(ds.samples[dist(gen)]);
  }
  return ds;
}


std::vector<double> evaluateLogCMF(FragmentLengthDistribution* fld) {
  std::vector<double> logPMFTemp;
  size_t minVal;
  size_t maxVal;
  fld->dumpPMF(logPMFTemp, minVal, maxVal);
  std::vector<double> logPMF(maxVal + 1, salmon::math::LOG_EPSILON);
  double sum = salmon::math::LOG_0;
  for (size_t i = 0; i < minVal; ++i) {
    sum = salmon::math::logAdd(sum, logPMF[i]);
  }
  for (size_t i = minVal; i < maxVal; ++i) {
    sum  = salmon::math::logAdd(sum, logPMFTemp[i-minVal]);
  }
  return fld->cmf(logPMF);
}


  LogCMFCache::LogCMFCache(FragmentLengthDistribution* fld,
                           bool singleEndLib, size_t fragUpdateThresh) :
    fld_(fld),
    singleEndLib_(singleEndLib),
    fragUpdateThresh_(fragUpdateThresh){
  }

  void LogCMFCache::refresh(size_t processedReads, bool burnedIn) {
      // The number of fragments we've seen since last time we re-cached the CMF
      size_t numNewFragments = static_cast<size_t>(processedReads - prevProcessedReads_);
      // We need to re-cache the CMF if this is our first pass through or if we've seen
      // at least fragUpdateThresh new fragments
      bool needFreshCMF = ((numNewFragments >= fragUpdateThresh_) or (prevProcessedReads_ == 0));
      // If we have a single-end library or we are burned in, then we don't need to
      // re-cache the CMF any more
      bool updateCachedCMF = (!singleEndLib_ and !burnedIn and needFreshCMF);

      // If we are going to attempt to model single mappings (part of a fragment)
      // then cache the FLD cumulative distribution for this mini-batch.  If
      // we are burned in or this is a single-end library, then the CMF is already
      // cached, so we don't need to worry about doing this work ourselves.
      if (updateCachedCMF) {
        cachedCMF_ = distribution_utils::evaluateLogCMF(fld_);
      }

      if ((numNewFragments >= fragUpdateThresh_) or (prevProcessedReads_ == 0)) {
        prevProcessedReads_ = processedReads;
      }
  }

  double LogCMFCache::getAmbigFragLengthProb(bool fwd, int32_t pos, int32_t rlen, int32_t tlen, bool burnedIn) const {
    // if this is forward, then the "virtual" mate would be downstream
    // if thsi is rc, then the "virtual" mate would be upstream
    int32_t maxFragLen = 0;
    int32_t sTxpLen = static_cast<int32_t>(tlen);
    if (fwd) {
      int32_t p1 = (pos < 0) ? 0 : pos;
      p1 = (p1 > sTxpLen) ? sTxpLen : p1;
      maxFragLen = (sTxpLen - p1);
    } else {
      int32_t p1 = pos + rlen;
      p1 = (p1 < 0) ? 0 : p1;
      p1 = (p1 > sTxpLen) ? sTxpLen : p1;
      maxFragLen = (p1);
    }

    bool useFLD = (singleEndLib_ or burnedIn);
    double refLengthCM = cmfValue_(static_cast<size_t>(tlen), useFLD);
    bool computeMass = !salmon::math::isLog0(refLengthCM);
    double maxLenProb = cmfValue_(maxFragLen, useFLD);
    return  (computeMass) ? (maxLenProb - refLengthCM) : salmon::math::LOG_EPSILON;
  }

  inline double LogCMFCache::cmfValue_(size_t len, bool useFLD) const {
    return (useFLD) ? fld_->cmf(len) :
      (len < cachedCMF_.size()) ? cachedCMF_[len] : cachedCMF_.back();
  }

} // namespace distribution_utils
