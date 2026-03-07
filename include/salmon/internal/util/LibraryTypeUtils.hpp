#ifndef SALMON_INTERNAL_UTIL_LIBRARY_TYPE_UTILS_HPP
#define SALMON_INTERNAL_UTIL_LIBRARY_TYPE_UTILS_HPP

#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "salmon/internal/model/LibraryFormat.hpp"
#include "salmon/internal/quant/ReadLibrary.hpp"

namespace salmon::utils {

LibraryFormat parseLibraryFormatStringNew(std::string& fmt);

std::vector<ReadLibrary>
extractReadLibraries(boost::program_options::parsed_options& orderedOptions);

LibraryFormat parseLibraryFormatString(std::string& fmt);

} // namespace salmon::utils

#endif
