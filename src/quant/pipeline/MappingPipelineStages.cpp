#include "salmon/internal/quant/pipeline/MappingPipelineStages.hpp"

#include <functional>

#include <Eigen/Dense>

#include "salmon/internal/inference/CollapsedEMOptimizer.hpp"
#include "salmon/internal/inference/CollapsedGibbsSampler.hpp"
#include "salmon/internal/model/Transcript.hpp"
#include "salmon/internal/model/TranscriptGroup.hpp"
#include "salmon/internal/util/LibraryTypeUtils.hpp"
#include "salmon/internal/util/QuantOptionsUtils.hpp"
#include "salmon/internal/util/SalmonUtils.hpp"

namespace salmon::pipeline::stages {

bool stageProcessMappingOptions(QuantPipelineContext& ctx) {
  ctx.sopt.quantMode = SalmonQuantMode::MAP;
  bool optionsOK =
      salmon::utils::processQuantOptions(ctx.sopt, ctx.vm, ctx.numBiasSamples);
  if (!optionsOK) {
    return false;
  }

  ctx.fileLog = ctx.sopt.fileLog;
  ctx.jointLog = ctx.sopt.jointLog;
  ctx.indexDirectory = ctx.sopt.indexDirectory;
  ctx.outputDirectory = ctx.sopt.outputDirectory;
  return true;
}

std::vector<ReadLibrary> stageParseReadLibraries(QuantPipelineContext& ctx) {
  ctx.jointLog->info("parsing read library format");
  return salmon::utils::extractReadLibraries(ctx.orderedOptions);
}

int stageFinalizeMappingOutputs(QuantPipelineContext& ctx,
                                MappingReadExperiment& experiment,
                                MappingStatistics& mstats, GZipWriter& gzw,
                                const boost::filesystem::path& outputDirectory) {
  namespace bfs = boost::filesystem;
  auto& sopt = ctx.sopt;
  auto jointLog = ctx.jointLog;

  if (!sopt.skipQuant) {
    CollapsedEMOptimizer optimizer;
    jointLog->info("Starting optimizer");
    salmon::utils::normalizeAlphas(sopt, experiment);
    bool optSuccess = optimizer.optimize(experiment, sopt, 0.01, 10000);

    if (!optSuccess) {
      jointLog->error(
          "The optimization algorithm failed. This is likely the result of "
          "bad input (or a bug). If you cannot track down the cause, please "
          "report this issue on GitHub.");
      return 1;
    }
    jointLog->info("Finished optimizer");

    jointLog->info("writing output \n");
    gzw.writeAbundances(sopt, experiment, false);

    if (sopt.numGibbsSamples > 0) {
      jointLog->info("Starting Gibbs Sampler");
      CollapsedGibbsSampler sampler;
      gzw.setSamplingPath(sopt);
      std::function<bool(const std::vector<double>&)> bsWriter =
          [&gzw](const std::vector<double>& alphas) -> bool {
        return gzw.writeBootstrap(alphas, true);
      };

      bool sampleSuccess =
          sampler.sample(experiment, sopt, bsWriter, sopt.numGibbsSamples);
      if (!sampleSuccess) {
        jointLog->error("Encountered error during Gibbs sampling.\n"
                        "This should not happen.\n"
                        "Please file a bug report on GitHub.\n");
        return 1;
      }
      jointLog->info("Finished Gibbs Sampler");
    } else if (sopt.numBootstraps > 0) {
      gzw.setSamplingPath(sopt);
      std::function<bool(const std::vector<double>&)> bsWriter =
          [&gzw](const std::vector<double>& alphas) -> bool {
        return gzw.writeBootstrap(alphas);
      };

      jointLog->info("Starting Bootstrapping");
      bool bootstrapSuccess =
          optimizer.gatherBootstraps(experiment, sopt, bsWriter, 0.01, 10000);
      jointLog->info("Finished Bootstrapping");
      if (!bootstrapSuccess) {
        jointLog->error("Encountered error during bootstrapping.\n"
                        "This should not happen.\n"
                        "Please file a bug report on GitHub.\n");
        return 1;
      }
    }

    if (ctx.vm.count("geneMap")) {
      try {
        auto geneOutputDirectory = outputDirectory;
        salmon::utils::generateGeneLevelEstimates(sopt.geneMapPath,
                                                  geneOutputDirectory);
      } catch (std::invalid_argument& e) {
        fmt::print(stderr,
                   "Error: [{}] when trying to compute gene-level "
                   "estimates. The gene-level file(s) may not exist",
                   e.what());
      }
    }
  } else if (sopt.dumpEqWeights) {
    jointLog->info("Finalizing combined weights for equivalence classes.");
    auto& eqVec = experiment.equivalenceClassBuilder().eqVec();
    bool noRichEq = sopt.noRichEqClasses;
    bool useEffectiveLengths = !sopt.noEffectiveLengthCorrection;
    std::vector<Transcript>& transcripts = experiment.transcripts();
    Eigen::VectorXd effLens(transcripts.size());

    for (size_t i = 0; i < transcripts.size(); ++i) {
      auto& txp = transcripts[i];
      effLens(i) = useEffectiveLengths
                       ? std::exp(txp.getCachedLogEffectiveLength())
                       : txp.RefLength;
    }

    for (size_t eqID = 0; eqID < eqVec.size(); ++eqID) {
      auto& kv = eqVec[eqID];
      const TranscriptGroup& k = kv.first;
      size_t classSize = kv.second.weights.size();
      auto& v = kv.second;
      double wsum{0.0};

      for (size_t i = 0; i < classSize; ++i) {
        auto tid = k.txps[i];
        double el = effLens(tid);
        if (el <= 1.0) {
          el = 1.0;
        }
        if (noRichEq) {
          v.weights[i] = 1.0;
        }
        auto probStartPos = 1.0 / el;
        double wt = sopt.eqClassMode ? v.weights[i]
                                     : v.count * v.weights[i] * probStartPos;
        v.combinedWeights.push_back(wt);
        wsum += wt;
      }

      double wnorm = 1.0 / wsum;
      for (size_t i = 0; i < classSize; ++i) {
        v.combinedWeights[i] = v.combinedWeights[i] * wnorm;
      }
    }
    jointLog->info("done.");
  }

  if (sopt.dumpEq) {
    jointLog->info("writing equivalence class counts.");
    gzw.writeEquivCounts(sopt, experiment);
    jointLog->info("done writing equivalence class counts.");
  }

  bfs::path libCountFilePath = outputDirectory / "lib_format_counts.json";
  experiment.summarizeLibraryTypeCounts(libCountFilePath);

  if (!sopt.noFragLengthDist) {
    bfs::path distFileName = sopt.paramsDirectory / "flenDist.txt";
    std::unique_ptr<std::FILE, int (*)(std::FILE*)> distOut(
        std::fopen(distFileName.c_str(), "w"), std::fclose);
    fmt::print(distOut.get(), "{}\n",
               experiment.fragmentLengthDistribution()->toString());
  }

  if (sopt.writeUnmappedNames) {
    auto l = sopt.unmappedLog.get();
    if (l) {
      l->flush();
      if (sopt.unmappedFile) {
        sopt.unmappedFile->close();
      }
    }
  }

  if (sopt.writeOrphanLinks) {
    auto l = sopt.orphanLinkLog.get();
    if (l) {
      l->flush();
      if (sopt.orphanLinkFile) {
        sopt.orphanLinkFile->close();
      }
    }
  }

  if (sopt.qmFileName != "") {
    sopt.qmLog->flush();
    if (sopt.qmFileName != "-") {
      sopt.qmFile.close();
    }
  }

  sopt.runStopTime = salmon::utils::getCurrentTimeAsString();
  gzw.writeMeta(sopt, experiment, mstats);
  spdlog::drop_all();
  return 0;
}

void stageDispatchByIndexType(QuantPipelineContext& ctx,
                              MappingReadExperiment& experiment,
                              std::vector<ReadLibrary>& readLibraries,
                              const std::function<void()>& runPuffQuant) {
  auto& sopt = ctx.sopt;
  auto jointLog = ctx.jointLog;
  auto indexType = experiment.getIndex()->indexType();

  switch (indexType) {
  case SalmonIndexType::FMD: {
    fmt::MemoryWriter infostr;
    infostr << "This version of salmon does not support FMD indexing.";
    throw std::invalid_argument(infostr.str());
  } break;
  case SalmonIndexType::QUASI: {
    fmt::MemoryWriter infostr;
    infostr << "This version of salmon does not support RapMap-based indexing.";
    throw std::invalid_argument(infostr.str());
  } break;
  case SalmonIndexType::PUFF: {
    if (sopt.gcBiasCorrect) {
      for (auto& rl : readLibraries) {
        if (rl.format().type != ReadType::PAIRED_END) {
          jointLog->warn(
              "Fragment GC bias correction is currently *experimental* "
              "in single-end libraries.  Please use this option "
              "with caution.");
        }
      }
    }
    runPuffQuant();
  } break;
  }
}

} // namespace salmon::pipeline::stages
