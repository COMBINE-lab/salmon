#include <cstdio>
#include <iostream>
#include <map>
#include <tuple>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/geometric.hpp>


#include "ONTAlignmentModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"

ONTAlignmentModel::ONTAlignmentModel(double alpha, uint32_t readBins)
    : isEnabled_(true)
    , readBins_(readBins)
    , errorModel_(maxReadLen / binLen + 1)
{
  for(auto& avg: errorModel_) {
    avg.number = 0;
    avg.sum    = 0.0;
  }
}

double ONTAlignmentModel::logLikelihood(const UnpairedRead& hit, const UnpairedRead& primary, Transcript& ref) {
  ErrorCount counts;
  if(!computeErrorCount(hit.read, primary.read, ref, counts, "logLikelihood")) {
    if(logger_)
      logger_->warn("in logLikelihood() error parsing CIGAR string");
    return salmon::math::LOG_1;
  }

  const uint32_t readLen   = alnLen(hit, primary);
  const uint32_t alignLen  = readLen - counts.clips;
  const double   errorRate = (double)counts.ims() / alignLen;

  int32_t bin = std::min(alignLen / binLen, (uint32_t)errorModel_.size() - 1);
  auto& average = errorModel_[bin];
  if(average.number == 0) { // Falls into a bin that is empty after training
    if(logger_)
      logger_->warn("read {} (length {} - align length {}) has no trained error model",
                    bam_name(hit.read), readLen, alignLen);
    return salmon::math::LOG_1;
  }

  // Binomial / Geometric distribution probability of an event (an
  // error in the sequence, whether a mismatch of indel).
  const double p = average.sum / average.number;

  // Likelihood, based on average number of errors to get a read
  // further away from mode (as number of errors).
  using boost::math::binomial;
  if(p <= 0) {
    if(logger_)
      logger_->warn("Negative probability!!!, read {}", bam_name(hit.read));
    return salmon::math::LOG_1;
  }

  binomial errorDist(alignLen, p);
  const int32_t errorMedian = median(errorDist);
  const int32_t offMedian = std::abs(errorMedian - counts.ims());
  const double errorLikelihood = cdf(errorDist, std::max(errorMedian - offMedian, (int32_t)0)) +
    cdf(complement(errorDist, std::min(errorMedian + offMedian, (int32_t)alignLen)));
  if(errorLikelihood <= 0.0) {
    if(logger_)
      logger_->warn("negative likelihood!!!, read {}", bam_name(hit.read));
  }

  // Likelihood to have so many bases soft clipped based on the
  // average error rate.
  // using boost::math::geometric;
  // geometric clipDist(1.0 - p);
  // const double clipLikelihood = 1.0 - cdf(clipDist, counts.clips);

  return log(errorLikelihood); // + log(clipLikelihood);
}

// Update the probability model. The reads are binned based on their
// length (the length after soft clipping). For each bin, we assume a
// Binomial distribution B(p,n) (where n is the length of the read).
//
// The probability p in each bin is estimated as the empirical mean of
// the number of errors in the reads (all reads counted the same at
// this point: indels and mutations).
void ONTAlignmentModel::update(const UnpairedRead& hit, const UnpairedRead& primary,
                               Transcript& ref, double p, double mass) {
  if (mass == salmon::math::LOG_0) {
    return;
  }
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return;
  }

  ErrorCount counts;
  if(!computeErrorCount(hit.read, primary.read, ref, counts, "update")) {
    if(logger_)
      logger_->warn("in update() error parsing CIGAR string");
    return;
  }

  // Not taking p and mass into account. What's up with those?
  const int32_t readLen   = alnLen(hit, primary);
  const int32_t alignLen  = readLen - counts.clips;
  const double  errorRate = (double)counts.ims() / alignLen;
  if(errorRate > 1.0) { // Should not happen
    if (logger_) {
      logger_->warn("(in update()) CIGAR string for read [{}] "
                    "seems inconsistent. It implied an error rate "
                    "greater than 1",
                    bam_name(hit.read));
    }
    return;
  }

  int32_t bin = std::min(alignLen / binLen, (uint32_t)errorModel_.size() - 1);
  ++errorModel_[bin].number;
  salmon::utils::incLoop(errorModel_[bin].sum, errorRate);
}

void ONTAlignmentModel::printModel(std::ostream& os) {
  for(size_t i = 0; i < errorModel_.size(); ++i) {
    if(errorModel_[i].number == 0) continue;
    const auto p = (errorModel_[i].sum / errorModel_[i].number);
    const auto n = i * binLen;
    os << (i * binLen) << " - " << ((i+1) * binLen) << ' ' << p
       << ' ' << errorModel_[i].number <<'\n';
  }
}
