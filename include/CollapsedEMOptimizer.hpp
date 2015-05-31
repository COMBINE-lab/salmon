#ifndef COLLAPSED_EM_OPTIMIZER_HPP
#define COLLAPSED_EM_OPTIMIZER_HPP

#include <unordered_map>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"

#include "cuckoohash_map.hh"
#include "Eigen/Dense"

class CollapsedEMOptimizer {
    public:
        using VecType = std::vector<tbb::atomic<double>>;
        CollapsedEMOptimizer();

        template <typename ExpT>
        bool optimize(ExpT& readExp,
                      SalmonOpts& sopt,
                      double tolerance = 0.01,
                      uint32_t maxIter = 1000);
    private:
        void EMUpdate_(
        std::vector<std::pair<const TranscriptGroup, TGValue>>& eqVec,
        std::vector<Transcript>& transcripts,
        Eigen::VectorXd& effLens,
        const VecType& alphaIn,
        VecType& alphaOut);

};

#endif // COLLAPSED_EM_OPTIMIZER_HPP
