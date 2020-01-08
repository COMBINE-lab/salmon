#ifndef PROGRAM_OPTIONS_GENERATOR_HPP
#define PROGRAM_OPTIONS_GENERATOR_HPP

#include <boost/program_options.hpp>
#include "SalmonOpts.hpp"

namespace salmon {
namespace po = boost::program_options;
class ProgramOptionsGenerator{
public:
  po::options_description getMappingInputOptions(SalmonOpts& sopt);
  po::options_description getAlignmentInputOptions(SalmonOpts& sopt);

  po::options_description getBasicOptions(SalmonOpts& sopt);

  po::options_description getMappingSpecificOptions(SalmonOpts& sopt);
  po::options_description getAlignmentSpecificOptions(SalmonOpts& sopt);
  po::options_description getAlevinBasicOptions(SalmonOpts& sopt);
  po::options_description getAlevinDevsOptions();
  po::options_description getAdvancedOptions(int32_t& numBiasSamples, SalmonOpts& sopt);
  po::options_description getHiddenOptions(SalmonOpts& sopt);
  po::options_description getTestingOptions(SalmonOpts& sopt);
  po::options_description getDeprecatedOptions(SalmonOpts& /*sopt*/);
};

}

#endif // PROGRAM_OPTIONS_GENERATOR_HPP
