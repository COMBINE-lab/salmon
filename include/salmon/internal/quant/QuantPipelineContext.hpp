#ifndef SALMON_INTERNAL_QUANT_QUANT_PIPELINE_CONTEXT_HPP
#define SALMON_INTERNAL_QUANT_QUANT_PIPELINE_CONTEXT_HPP

#include <cstdint>
#include <string>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <spdlog/spdlog.h>

#include "salmon/internal/config/SalmonOpts.hpp"

namespace salmon::pipeline {

struct QuantPipelineContext {
  QuantPipelineContext(SalmonOpts& opts,
                       boost::program_options::variables_map& variables,
                       boost::program_options::parsed_options& ordered,
                       int32_t biasSamples)
      : sopt(opts),
        vm(variables),
        orderedOptions(ordered),
        numBiasSamples(biasSamples) {}

  SalmonOpts& sopt;
  boost::program_options::variables_map& vm;
  boost::program_options::parsed_options& orderedOptions;
  int32_t numBiasSamples{0};

  std::shared_ptr<spdlog::logger> jointLog;
  std::shared_ptr<spdlog::logger> fileLog;
  boost::filesystem::path indexDirectory;
  boost::filesystem::path outputDirectory;
  std::string commandComment;
};

} // namespace salmon::pipeline

#endif
