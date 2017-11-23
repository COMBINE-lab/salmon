#ifndef __SALMON_EXCEPTIONS_HPP__
#define __SALMON_EXCEPTIONS_HPP__

#include <sstream>
#include <stdexcept>

class InsufficientAssignedFragments : public std::runtime_error {
public:
  InsufficientAssignedFragments(uint64_t numAssigned,
                                uint64_t numRequired) noexcept;
  InsufficientAssignedFragments(InsufficientAssignedFragments&& other);
  virtual const char* what() const throw();

private:
  uint64_t numAssigned_;
  uint64_t numRequired_;
  std::ostringstream cnvt;
  std::string msg_;
};

#endif //__SALMON_EXCEPTIONS_HPP__
