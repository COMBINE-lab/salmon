#ifndef COLLAPSED_EM_OPTIMIZER_HPP
#define COLLAPSED_EM_OPTIMIZER_HPP

#include <functional>
#include <unordered_map>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"

#include "Eigen/Dense"
#include "cuckoohash_map.hh"

class BootstrapWriter;

class CollapsedEMOptimizer {
public:
  using VecType = std::vector<tbb::atomic<double>>;
  using SerialVecType = std::vector<double>;
  CollapsedEMOptimizer();

  template <typename ExpT>
  bool optimize(
      ExpT& readExp, SalmonOpts& sopt,
      double tolerance =
          0.01, // A EM termination criteria, adopted from Bray et al. 2016
      uint32_t maxIter =
          1000); // A EM termination criteria, adopted from Bray et al. 2016

  template <typename ExpT>
  bool gatherBootstraps(
      ExpT& readExp, SalmonOpts& sopt,
      std::function<bool(const std::vector<double>&)>& writeBootstrap,
      double relDiffTolerance, uint32_t maxIter);
};

#endif // COLLAPSED_EM_OPTIMIZER_HPP
