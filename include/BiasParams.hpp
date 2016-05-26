#ifndef __BIAS_PARAMS__
#define __BIAS_PARAMS__

#include "SBModel.hpp"
#include "ReadKmerDist.hpp"
#include "SalmonMath.hpp"
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
    std::vector<double> observedGCMass = std::vector<double>(101, salmon::math::LOG_0);

    ReadKmerDist<8, uint32_t> seqBiasFW;
    ReadKmerDist<8, uint32_t> seqBiasRC;

  /**
   * Sequence-specific bias models
   **/
    SBModel seqBiasModelFW;
    SBModel seqBiasModelRC;

    BiasParams(bool seqBiasPseudocount=false) : seqBiasFW(seqBiasPseudocount), seqBiasRC(seqBiasPseudocount) {}
};

#endif //__GC_BIAS_PARAMS__
