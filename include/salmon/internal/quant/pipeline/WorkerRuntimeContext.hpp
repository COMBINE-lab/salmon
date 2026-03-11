#ifndef SALMON_INTERNAL_QUANT_PIPELINE_WORKER_RUNTIME_CONTEXT_HPP
#define SALMON_INTERNAL_QUANT_PIPELINE_WORKER_RUNTIME_CONTEXT_HPP

#include <cmath>
#include <cstdint>
#include <vector>

#include "salmon/internal/quant/BiasParams.hpp"
#include "salmon/internal/util/SalmonMath.hpp"
#include "salmon/internal/util/SalmonUtils.hpp"

namespace salmon::pipeline {

struct WorkerRuntimeContext {
  WorkerRuntimeContext(uint32_t workerThreadCount, size_t numConditionalGCBins,
                       size_t numFragGCBins)
      : workerThreads(workerThreadCount),
        observedBias(workerThreadCount,
                     BiasParams(numConditionalGCBins, numFragGCBins, false)) {}

  template <typename QuantContextT>
  void mergeObservedBias(QuantContextT& quantContext,
                         bool finalizePositionalBias) const {
    double gcFracFwd{0.0};
    double globalMass{salmon::math::LOG_0};
    double globalFwdMass{salmon::math::LOG_0};
    auto& globalGCMass = quantContext.observedGC();

    for (auto& localBias : observedBias) {
      globalGCMass.combineCounts(localBias.observedGCMass);

      auto& fw =
          quantContext.readBiasModelObserved(salmon::utils::Direction::FORWARD);
      auto& rc = quantContext.readBiasModelObserved(
          salmon::utils::Direction::REVERSE_COMPLEMENT);

      fw.combineCounts(localBias.seqBiasModelFW);
      rc.combineCounts(localBias.seqBiasModelRC);

      auto& posBiasesFW = quantContext.posBias(salmon::utils::Direction::FORWARD);
      auto& posBiasesRC =
          quantContext.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      for (size_t i = 0; i < posBiasesFW.size(); ++i) {
        posBiasesFW[i].combine(localBias.posBiasFW[i]);
        posBiasesRC[i].combine(localBias.posBiasRC[i]);
      }

      globalMass = salmon::math::logAdd(globalMass, localBias.massFwd);
      globalMass = salmon::math::logAdd(globalMass, localBias.massRC);
      globalFwdMass = salmon::math::logAdd(globalFwdMass, localBias.massFwd);
    }

    globalGCMass.normalize();
    if (globalMass != salmon::math::LOG_0) {
      if (globalFwdMass != salmon::math::LOG_0) {
        gcFracFwd = std::exp(globalFwdMass - globalMass);
      }
      quantContext.setGCFracForward(gcFracFwd);
    }

    if (finalizePositionalBias) {
      auto& posBiasesFW = quantContext.posBias(salmon::utils::Direction::FORWARD);
      auto& posBiasesRC =
          quantContext.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      for (size_t i = 0; i < posBiasesFW.size(); ++i) {
        posBiasesFW[i].finalize();
        posBiasesRC[i].finalize();
      }
    }
  }

  uint32_t workerThreads{0};
  std::vector<BiasParams> observedBias;
};

} // namespace salmon::pipeline

#endif
