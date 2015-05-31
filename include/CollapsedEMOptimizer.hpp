#ifndef COLLAPSED_EM_OPTIMIZER_HPP
#define COLLAPSED_EM_OPTIMIZER_HPP

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"

#include <unordered_map>

class CollapsedEMOptimizer {
    public:
        CollapsedEMOptimizer();

        template <typename ExpT>
        bool optimize(ExpT& readExp,
                      SalmonOpts& sopt,
                      double tolerance = 0.01,
                      uint32_t maxIter = 1000);
};

#endif // COLLAPSED_EM_OPTIMIZER_HPP
