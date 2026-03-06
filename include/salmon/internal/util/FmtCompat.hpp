#ifndef SALMON_INTERNAL_UTIL_FMT_COMPAT_HPP
#define SALMON_INTERNAL_UTIL_FMT_COMPAT_HPP

#include <sstream>
#include <string>

namespace fmt {

// Compatibility shim for legacy fmt::MemoryWriter usage while migrating to
// newer bundled fmt versions shipped with modern spdlog.
class MemoryWriter {
public:
  MemoryWriter() = default;

  template <typename T>
  MemoryWriter& operator<<(const T& value) {
    stream_ << value;
    return *this;
  }

  MemoryWriter& operator<<(std::ostream& (*manip)(std::ostream&)) {
    stream_ << manip;
    return *this;
  }

  std::string str() const { return stream_.str(); }

  void clear() {
    stream_.str("");
    stream_.clear();
  }

private:
  std::ostringstream stream_;
};

} // namespace fmt

#endif
