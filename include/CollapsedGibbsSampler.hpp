#ifndef COLLAPSED_GIBBS_SAMPLER_HPP
#define COLLAPSED_GIBBS_SAMPLER_HPP

#include <unordered_map>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include "SalmonOpts.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

class BootstrapWriter;

class CollapsedGibbsSampler {
    public:
        using VecType = std::vector<double>;
        CollapsedGibbsSampler();

        template <typename ExpT>
        bool sample(ExpT& readExp,
                      SalmonOpts& sopt,
                      BootstrapWriter* bwriter,
                      uint32_t numSamples = 500);

};

#endif // COLLAPSED_EM_OPTIMIZER_HPP

