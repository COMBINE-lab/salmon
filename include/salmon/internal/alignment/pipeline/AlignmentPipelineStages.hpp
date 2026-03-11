#ifndef SALMON_INTERNAL_ALIGNMENT_PIPELINE_ALIGNMENT_PIPELINE_STAGES_HPP
#define SALMON_INTERNAL_ALIGNMENT_PIPELINE_ALIGNMENT_PIPELINE_STAGES_HPP

#include <functional>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "salmon/internal/model/LibraryFormat.hpp"
#include "salmon/internal/quant/QuantPipelineContext.hpp"

namespace salmon::pipeline::stages {

bool stageProcessAlignmentOptions(QuantPipelineContext& ctx);

std::vector<boost::filesystem::path> stageCollectAlignmentInputs(
    const std::vector<std::string>& alignmentFileNames,
    const std::string& eqclassesFileName, bool hasAlignments, bool hasEqclasses);

int stageFinalizeAlignmentOutputs(QuantPipelineContext& ctx, bool success,
                                  const boost::filesystem::path& outputDirectory);

bool stageDispatchByLibraryType(
    QuantPipelineContext& ctx, LibraryFormat& libFmt, bool hasEqclasses,
    const std::function<bool()>& runSingleEndEqClasses,
    const std::function<bool()>& runSingleEnd,
    const std::function<bool()>& runSingleEndONT,
    const std::function<bool()>& runPairedEnd);

} // namespace salmon::pipeline::stages

#endif
