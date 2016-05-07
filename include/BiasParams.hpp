#ifndef __BIAS_PARAMS__
#define __BIAS_PARAMS__

#include "SBModel.hpp"
#include "ReadKmerDist.hpp"
#include "SalmonMath.hpp"
#include <vector>

struct BiasParams {
    double massFwd{salmon::math::LOG_0};
    double massRC{salmon::math::LOG_0};
    std::vector<double> observedGCMass = std::vector<double>(101, salmon::math::LOG_0);
    ReadKmerDist<8, uint32_t> seqBiasFW;
    ReadKmerDist<8, uint32_t> seqBiasRC;

    SBModel seqBiasModelFW;
    SBModel seqBiasModelRC;

    BiasParams(bool seqBiasPseudocount=false) : seqBiasFW(seqBiasPseudocount), seqBiasRC(seqBiasPseudocount) {}
};

#endif //__GC_BIAS_PARAMS__
