#include "SalmonExceptions.hpp"


InsufficientAssignedFragments::InsufficientAssignedFragments(uint64_t numAssigned, uint64_t numRequired) noexcept :
  std::runtime_error(""), numAssigned_(numAssigned), numRequired_(numRequired) {
  cnvt.str("");
  cnvt << "Salmon was only able to assign " << numAssigned_ << " fragments to transcripts in the index, "
    "but the minimum number of required assigned fragments (--minAssignedFrags) "
    "was " << numRequired_ << ". This could be indicative of a mismatch between the reference and "
    "sample, or a very bad sample.  You can change the --minAssignedFrags "
    "parameter to force Salmon to quantify with fewer assigned fragments (must have at least 1).";
  msg_ = cnvt.str();
}

const char* InsufficientAssignedFragments::what() const throw() {
  return msg_.c_str();
}
