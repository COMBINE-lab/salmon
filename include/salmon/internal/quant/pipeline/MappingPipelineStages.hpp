#ifndef SALMON_INTERNAL_QUANT_PIPELINE_MAPPING_PIPELINE_STAGES_HPP
#define SALMON_INTERNAL_QUANT_PIPELINE_MAPPING_PIPELINE_STAGES_HPP

#include <functional>
#include <vector>

#include <boost/filesystem.hpp>

#include "salmon/internal/model/Transcript.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "salmon/internal/output/GZipWriter.hpp"
#include "salmon/internal/output/MappingStatistics.hpp"
#include "salmon/internal/quant/QuantPipelineContext.hpp"
#include "salmon/internal/quant/ReadExperiment.hpp"
#include "salmon/internal/quant/ReadLibrary.hpp"

namespace salmon::pipeline::stages {

using MappingReadExperiment = ReadExperiment<EquivalenceClassBuilder<TGValue>>;

bool stageProcessMappingOptions(QuantPipelineContext& ctx);

std::vector<ReadLibrary> stageParseReadLibraries(QuantPipelineContext& ctx);

int stageFinalizeMappingOutputs(QuantPipelineContext& ctx,
                                MappingReadExperiment& experiment,
                                MappingStatistics& mstats, GZipWriter& gzw,
                                const boost::filesystem::path& outputDirectory);

void stageDispatchByIndexType(QuantPipelineContext& ctx,
                              MappingReadExperiment& experiment,
                              std::vector<ReadLibrary>& readLibraries,
                              const std::function<void()>& runPuffQuant);

} // namespace salmon::pipeline::stages

#endif
