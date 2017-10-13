#ifndef __BIAS_PARAMS__
#define __BIAS_PARAMS__

#include "DistributionUtils.hpp"
#include "GCFragModel.hpp"
#include "ReadKmerDist.hpp"
#include "SBModel.hpp"
#include "SalmonMath.hpp"
#include "SimplePosBias.hpp"
#include <vector>

struct BiasParams {
  double massFwd{salmon::math::LOG_0};
  double massRC{salmon::math::LOG_0};

  /**
   * Positional bias
   **/
  std::vector<SimplePosBias> posBiasFW;
  std::vector<SimplePosBias> posBiasRC;

  /**
   * fragment-GC bias counts
   **/
  // std::vector<double> observedGCMass = std::vector<double>(101,
  // salmon::math::LOG_0);
  GCFragModel observedGCMass;

  ReadKmerDist<8, uint32_t> seqBiasFW;
  ReadKmerDist<8, uint32_t> seqBiasRC;

  /**
   * Sequence-specific bias models
   **/
  SBModel seqBiasModelFW;
  SBModel seqBiasModelRC;

  BiasParams(size_t numCondBins = 3, size_t numGCBins = 101,
             bool seqBiasPseudocount = false)
      : seqBiasFW(seqBiasPseudocount), seqBiasRC(seqBiasPseudocount),
        posBiasFW(5), posBiasRC(5), observedGCMass(numCondBins, numGCBins) {}
};

#endif //__GC_BIAS_PARAMS__
