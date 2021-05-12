#include <cstdio>
#include <iostream>
#include <map>
#include <tuple>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/geometric.hpp>

#include "dbg.hpp"
#include "ONTAlignmentModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"

ONTAlignmentModel::ONTAlignmentModel(double alpha, uint32_t readBins)
    : isEnabled_(true)
    , readBins_(readBins)
    , printed(false)
    , errorModel_(maxReadLen / binLen + 1)
    , frontClipModel_(maxReadLen / binLen + 1)
    , backClipModel_(maxReadLen / binLen + 1)
{ }

double ONTAlignmentModel::logLikelihood(const UnpairedRead& hit, const UnpairedRead& primary, Transcript& ref) {
  using salmon::math::LOG_0;
  using salmon::math::LOG_1;
  constexpr auto dmin = std::numeric_limits<double>::min();
  constexpr double llMin = 1e-10; // Ignore alignment with likelihood below that number

  // if(!printed) {
  //   std::unique_lock<std::mutex> lock(outputMutex_, std::try_to_lock);
  //   if(lock.owns_lock()) {
  //       printed = true;
  //       printModel(std::cerr);
  //     }
  // }

  // If chimeric alignment, doesn't align to this transcript. Return 0
  // probability.
  if(bam_aux_find(hit.read, "SA"))
    return LOG_0;

  ErrorCount counts;
  if(!computeErrorCount(hit.read, primary.read, ref, counts, "logLikelihood")) {
    if(logger_)
      logger_->warn("in logLikelihood() error parsing CIGAR string");
    return LOG_1;
  }

  const uint32_t cigarRLen = alnLen(hit, primary); // Read length minus hard clips
  if(counts.sclips() >= cigarRLen)
    return LOG_0; // Empty alignment!

  const uint32_t alignLen  = cigarRLen - counts.clips(); // Length of aligned part (no soft clip)
  const double   errorRate = (double)counts.ims() / alignLen;
  const int32_t  errorBin  = std::min(alignLen / binLen, (uint32_t)errorModel_.size() - 1);
  const int32_t  frontClipBin   = std::min(cigarRLen / binLen, (uint32_t)frontClipModel_.size() - 1);
  const int32_t  backClipBin    = std::min(cigarRLen / binLen, (uint32_t)backClipModel_.size() - 1);
  const auto&    errorAvg  = errorModel_[errorBin];
  const auto&    frontClipAvg   = frontClipModel_[frontClipBin];
  const auto&    backClipAvg    = backClipModel_[backClipBin];

  double errorllh = LOG_1, frontClipllh = LOG_1, backClipllh = LOG_1; // Error and clip Log Likelihood (front and back)

  if(errorAvg.mass > dmin && frontClipAvg.mass > dmin && backClipAvg.mass > dmin) {
    // Binomial / Geometric distribution probability of an event (an
    // error in the sequence, whether a mismatch of indel).
    const double errorP = errorAvg.sum / errorAvg.mass;

    // Likelihood, based on average mass of errors to get a read
    // further away from mode (as number of errors).
    using boost::math::binomial;
    binomial      errorDist(alignLen, errorP);
    const int32_t errorMedian     = median(errorDist);
    const int32_t offMedian       = std::abs(errorMedian - counts.ims());
    const double  errorLikelihood =  cdf(errorDist, std::max(errorMedian - offMedian, (int32_t)0)) +
      1.0 - cdf(complement(errorDist, std::min(errorMedian + offMedian, (int32_t)alignLen)));
    errorllh = errorLikelihood < llMin ? LOG_0 : std::log(errorLikelihood);
  } else { // Falls into a bin that is empty after training
    if(logger_)
      logger_->warn("read {} (length {} - align length {}) has no trained error model",
                    bam_name(hit.read), cigarRLen, alignLen);
  }

  // Likelihood to have so many bases soft clipped based on the
  // average error rate. Don't penalize for having fewer clipped bases
  // than average, only if more.
  // front clips:
  if(frontClipAvg.sum > dmin && frontClipAvg.mass > dmin) {
    using        boost::math::geometric;
    const double  mean          = frontClipAvg.sum / frontClipAvg.mass;
    geometric     clipDist(1.0 / (mean + 1.0));
    const int32_t rmean         = std::round(mean);
    const auto    clips         = counts.fclips();
    const double clipLikelihood =
      clips <= rmean
      ? 1.0
      : (1.0 - cdf(clipDist, clips)) / (1.0 - cdf(clipDist, rmean));
    frontClipllh = clipLikelihood < llMin ? LOG_0 : std::log(clipLikelihood);
  } else {
    if(logger_)
      logger_->warn("read {} (length {}) has no trained clipping model",
                    bam_name(hit.read), cigarRLen);
  }

  // back clips:
  if(backClipAvg.sum > dmin && backClipAvg.mass > dmin) {
    using        boost::math::geometric;
    const double  mean          = backClipAvg.sum / backClipAvg.mass;
    geometric     clipDist(1.0 / (mean + 1.0));
    const int32_t rmean         = std::round(mean);
    const auto    clips         = counts.bclips();
    const double clipLikelihood =
      clips <= rmean
      ? 1.0
      : (1.0 - cdf(clipDist, clips)) / (1.0 - cdf(clipDist, rmean));
    backClipllh = clipLikelihood < llMin ? LOG_0 : std::log(clipLikelihood);
  } else {
    if(logger_)
      logger_->warn("read {} (length {}) has no trained clipping model",
                    bam_name(hit.read), cigarRLen);
  }

  return errorllh + frontClipllh + backClipllh;
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

  if(bam_aux_find(hit.read, "SA")) // Chimeric alignment. Ignore
    return;

  ErrorCount counts;
  if(!computeErrorCount(hit.read, primary.read, ref, counts, "update")) {
    if(logger_)
      logger_->warn("in update() error parsing CIGAR string");
    return;
  }

  // Update error model
  // Not taking p and mass into account. What's up with those?
  const int32_t readLen   = alnLen(hit, primary);
  const int32_t alignLen  = readLen - counts.sclips();
  const double  errorRate = (double)counts.ims() / alignLen;
  const double  clipRateFront  = (double)counts.fclips() / (readLen + counts.hclips());
  const double  clipRateBack   = (double)counts.bclips() / (readLen + counts.hclips());
  if(errorRate > 1.0 || clipRateFront > 1.0 || clipRateBack > 1.0) { // Should not happen
    if (logger_) {
      logger_->warn("(in update()) CIGAR string for read [{}] "
                    "seems inconsistent. It implied an error rate "
                    "greater than 1: {} {} {}",
                    bam_name(hit.read), errorRate, clipRateFront, clipRateBack);
    }
    return;
  }

  const double newMass = mass;
  { int32_t binIndex = std::min(alignLen / binLen, (uint32_t)errorModel_.size() - 1);
    auto& bin = errorModel_[binIndex];
    salmon::utils::incLoop(bin.mass, newMass);
    salmon::utils::incLoop(bin.sum, newMass * errorRate);
  }

  // Update front clip model
  { int32_t binIndex = std::min(readLen / binLen, (uint32_t)frontClipModel_.size() - 1);
    auto& bin = frontClipModel_[binIndex];
    salmon::utils::incLoop(bin.mass, newMass);
    salmon::utils::incLoop(bin.sum, (binIndex + 1) * binLen * newMass * clipRateFront);
  }

  // Update back clip model
  { int32_t binIndex = std::min(readLen / binLen, (uint32_t)backClipModel_.size() - 1);
    auto& bin = backClipModel_[binIndex];
    salmon::utils::incLoop(bin.mass, newMass);
    salmon::utils::incLoop(bin.sum, (binIndex + 1) * binLen * newMass * clipRateBack);
  }
}

void ONTAlignmentModel::printModel(std::ostream& os) {
  // dbg d(os);

  // d << "Model\n";
  // for(size_t i = 0; i < std::max(errorModel_.size(), frontClipModel_.size()); ++i) {
  //   const auto errorP = errorModel_[i].mass != 0.0 ? (errorModel_[i].sum / errorModel_[i].mass) : 0.0;
  //   const auto fclipP = frontClipModel_[i].mass != 0.0 ? (frontClipModel_[i].sum / frontClipModel_[i].mass) : 0.0;
  //   const auto bclipP = backClipModel_[i].mass != 0.0 ? (backClipModel_[i].sum / backClipModel_[i].mass) : 0.0;
  //   if(errorP == 0.0 && fclipP == 0.0 && bclipP == 0.0) continue;

  //   const auto n = i * binLen;
  //   d << (i * binLen) << " - " << ((i+1) * binLen) 
  //     << ' ' << errorP << ' ' << errorModel_[i].mass.load()
  //     << ' ' << fclipP << ' ' << frontClipModel_[i].mass.load()
  //     << ' ' << bclipP << ' ' << backClipModel_[i].mass.load() << '\n';
  // }
  // d << "--------------\n";
}
