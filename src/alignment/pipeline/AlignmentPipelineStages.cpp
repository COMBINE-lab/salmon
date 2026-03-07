#include "salmon/internal/alignment/pipeline/AlignmentPipelineStages.hpp"

#include <cstdlib>
#include <sstream>
#include <stdexcept>

#include <spdlog/fmt/fmt.h>

#include "salmon/internal/util/QuantOptionsUtils.hpp"
#include "salmon/internal/util/SalmonUtils.hpp"

namespace salmon::pipeline::stages {

bool stageProcessAlignmentOptions(QuantPipelineContext& ctx) {
  ctx.sopt.alnMode = true;
  ctx.sopt.quantMode = SalmonQuantMode::ALIGN;
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

std::vector<boost::filesystem::path> stageCollectAlignmentInputs(
    const std::vector<std::string>& alignmentFileNames,
    const std::string& eqclassesFileName, bool hasAlignments, bool hasEqclasses) {
  namespace bfs = boost::filesystem;
  std::vector<bfs::path> alignmentFiles;

  if (hasAlignments) {
    for (auto& alignmentFileName : alignmentFileNames) {
      bfs::path alignmentFile(alignmentFileName);
      if (!bfs::exists(alignmentFile)) {
        std::stringstream ss;
        ss << "The provided alignment file: " << alignmentFile
           << " does not exist!\n";
        throw std::invalid_argument(ss.str());
      }
      alignmentFiles.push_back(alignmentFile);
    }
  } else if (hasEqclasses) {
    bfs::path alignmentFile(eqclassesFileName);
    if (!bfs::exists(alignmentFile)) {
      std::stringstream ss;
      ss << "The provided eqclasses file: " << alignmentFile
         << " does not exist!\n";
      throw std::invalid_argument(ss.str());
    }
    alignmentFiles.push_back(alignmentFile);
  }

  return alignmentFiles;
}

int stageFinalizeAlignmentOutputs(QuantPipelineContext& ctx, bool success,
                                  const boost::filesystem::path& outputDirectory) {
  auto jointLog = ctx.jointLog;
  auto& sopt = ctx.sopt;

  if (!success) {
    jointLog->error(
        "Quantification was un-successful.  Please check the log "
        "for information about why quantification failed. If this "
        "problem persists, please report this issue on GitHub.");
    return 1;
  }

  if (ctx.vm.count("geneMap")) {
    try {
      auto geneOutputDirectory = outputDirectory;
      salmon::utils::generateGeneLevelEstimates(sopt.geneMapPath,
                                                geneOutputDirectory);
    } catch (std::exception& e) {
      fmt::print(stderr,
                 "Error: [{}] when trying to compute gene-level "
                 "estimates. The gene-level file(s) may not exist",
                 e.what());
    }
  }

  return 0;
}

bool stageDispatchByLibraryType(
    QuantPipelineContext& ctx, LibraryFormat& libFmt, bool hasEqclasses,
    const std::function<bool()>& runSingleEndEqClasses,
    const std::function<bool()>& runSingleEnd,
    const std::function<bool()>& runSingleEndONT,
    const std::function<bool()>& runPairedEnd) {
  auto& sopt = ctx.sopt;
  auto jointLog = ctx.jointLog;
  bool success{false};

  switch (libFmt.type) {
  case ReadType::SINGLE_END: {
    if (sopt.gcBiasCorrect) {
      jointLog->warn(
          "Fragment GC bias correction is currently *experimental* "
          "in single-end libraries.  Please use this option "
          "with caution.");
    }

    if (hasEqclasses) {
      success = runSingleEndEqClasses();
    } else {
      success = sopt.oxfordNanoporeModel ? runSingleEndONT() : runSingleEnd();
    }
  } break;
  case ReadType::PAIRED_END: {
    if (hasEqclasses) {
      jointLog->error(" Cannot quantify eqclasses in mode \n"
                      " Please report this on github");
      std::exit(1);
    }
    success = runPairedEnd();
  } break;
  default:
    std::stringstream errfmt;
    errfmt << "Cannot quantify library of unknown format " << libFmt;
    jointLog->error(errfmt.str());
    jointLog->flush();
    std::exit(1);
  }

  return success;
}

} // namespace salmon::pipeline::stages
