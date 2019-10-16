#ifndef COLLAPSED_GIBBS_SAMPLER_HPP
#define COLLAPSED_GIBBS_SAMPLER_HPP

#include <functional>
#include <unordered_map>

#include "tbb/atomic.h"
#include "tbb/task_scheduler_init.h"

#include "SalmonOpts.hpp"

#include "Eigen/Dense"
#include "cuckoohash_map.hh"

class BootstrapWriter;

class CollapsedGibbsSampler {
public:
  using VecType = std::vector<double>;
  CollapsedGibbsSampler();

  template <typename ExpT>
  bool sample(ExpT& readExp, SalmonOpts& sopt,
              std::function<bool(const std::vector<double>&)>& writeBootstrap,
              uint32_t numSamples = 500);

  template <typename ExpT>
  bool sampleAllele(
              ExpT& readExp, SalmonOpts& sopt,
              std::function<bool(const std::vector<double>&)>& writeBootstrap,
              std::unordered_map<std::string, size_t>& mergedTranscriptNameMap,
              std::unordered_map<size_t, size_t>& alleleToMergeTranacriptMapping,
              std::vector<Transcript>& mergedTranscripts,
              uint32_t numSamples = 500);
  /*
        template <typename ExpT>
        bool sampleMultipleChains(ExpT& readExp,
              SalmonOpts& sopt,
              std::function<bool(const std::vector<double>&)>& writeBootstrap,
              uint32_t numSamples = 500);
  */
};

#endif // COLLAPSED_EM_OPTIMIZER_HPP
