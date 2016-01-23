#ifndef __GC_BIAS_PARAMS__
#define __GC_BIAS_PARAMS__

#include "SalmonMath.hpp"
#include <vector>

struct GCBiasParams {
    double massFwd{salmon::math::LOG_0};
    double massRC{salmon::math::LOG_0};
    std::vector<double> observedGCMass = std::vector<double>(101, salmon::math::LOG_0);
};

#endif //__GC_BIAS_PARAMS__
