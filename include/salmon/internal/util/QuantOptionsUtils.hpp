#ifndef SALMON_INTERNAL_UTIL_QUANT_OPTIONS_UTILS_HPP
#define SALMON_INTERNAL_UTIL_QUANT_OPTIONS_UTILS_HPP

#include <cstdint>

#include <boost/program_options.hpp>

#include "SalmonOpts.hpp"

namespace salmon::utils {

bool validateOptionsAlignment_(SalmonOpts& sopt,
                               boost::program_options::variables_map& vm);

bool validateOptionsMapping_(SalmonOpts& sopt,
                             boost::program_options::variables_map& vm);

bool createAuxMapLoggers_(SalmonOpts& sopt,
                          boost::program_options::variables_map& vm);

bool processQuantOptions(SalmonOpts& sopt,
                         boost::program_options::variables_map& vm,
                         int32_t numBiasSamples);

} // namespace salmon::utils

#endif
