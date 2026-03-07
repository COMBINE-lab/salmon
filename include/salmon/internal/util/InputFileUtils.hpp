#ifndef SALMON_INTERNAL_UTIL_INPUT_FILE_UTILS_HPP
#define SALMON_INTERNAL_UTIL_INPUT_FILE_UTILS_HPP

#include <string>
#include <vector>

#include <boost/filesystem.hpp>

namespace salmon::utils {

bool peekBAMIsPaired(const boost::filesystem::path& fname);

size_t numberOfReadsInFastaFile(const std::string& fname);

bool readKmerOrder(const std::string& fname, std::vector<uint64_t>& kmers);

} // namespace salmon::utils

#endif
