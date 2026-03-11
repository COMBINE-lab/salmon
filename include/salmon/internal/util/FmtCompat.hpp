#ifndef SALMON_INTERNAL_UTIL_FMT_COMPAT_HPP
#define SALMON_INTERNAL_UTIL_FMT_COMPAT_HPP

#include <sstream>
#include <string>
#include <type_traits>

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

namespace salmon::fmtcompat {

template <typename Int,
          typename = std::enable_if_t<std::is_integral_v<Int>>>
inline std::string group_digits(Int value) {
  using UnsignedInt = std::make_unsigned_t<Int>;
  const bool negative = std::is_signed_v<Int> && (value < 0);
  const UnsignedInt magnitude = negative
                                    ? static_cast<UnsignedInt>(-(value + 1)) + 1
                                    : static_cast<UnsignedInt>(value);

  std::string raw = std::to_string(magnitude);
  for (std::ptrdiff_t i = static_cast<std::ptrdiff_t>(raw.size()) - 3; i > 0; i -= 3) {
    raw.insert(static_cast<std::size_t>(i), ",");
  }

  if (negative) {
    raw.insert(raw.begin(), '-');
  }
  return raw;
}

} // namespace salmon::fmtcompat

#endif
