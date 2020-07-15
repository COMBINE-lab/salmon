#include <cstdio>
#include <iostream>
#include <map>
#include <tuple>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY

#include "ONTAlignmentModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"

ONTAlignmentModel::ONTAlignmentModel(double alpha, uint32_t readBins)
    : isEnabled_(true), readBins_(readBins)
{ }

// ONTAlignmentModel::AlnModelProb ONTAlignmentModel::logLikelihood(bam_seq_t* read, Transcript& ref,
//                                                                  std::vector<AtomicMatrix<double>>& transitionProbs)
// {
//   throw std::runtime_error("ONTAlignmentModel not yet implemented");
//   return {0.0, 0.0};
// }

double ONTAlignmentModel::logLikelihood(const UnpairedRead& hit, Transcript& ref) {
  throw std::runtime_error("ONTAlignmentModel not yet implemented");
  double logLike = salmon::math::LOG_1;
  double bg = salmon::math::LOG_1;
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return logLike;
  }
  return logLike - bg;
}

void ONTAlignmentModel::update(const UnpairedRead& hit, Transcript& ref, double p,
                               double mass) {
  throw std::runtime_error("ONTAlignmentModel not yet implemented");
  if (mass == salmon::math::LOG_0) {
    return;
  }
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return;
  }
}

// void ONTAlignmentModel::update(
//     bam_seq_t* read, Transcript& ref, double p, double mass,
//     std::vector<AtomicMatrix<double>>& transitionProbs) {
//   throw std::runtime_error("ONTAlignmentModel not yet implemented");
// }

