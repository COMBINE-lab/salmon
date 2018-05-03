#ifndef SALMON_MATH_HPP
#define SALMON_MATH_HPP

#include "boost/config.hpp"

// If we have built-ins, do as Boost does

#ifndef BOOST_LIKELY
#if defined(__has_builtin)
#if __has_builtin(__builtin_expect)
#define BOOST_LIKELY(x) __builtin_expect(x, 1)
#endif
#endif
#endif

#ifndef BOOST_UNLIKELY
#if defined(__has_builtin)
#if __has_builtin(__builtin_expect)
#define BOOST_UNLIKELY(x) __builtin_expect(x, 0)
#endif
#endif
#endif

// If we didn't have those built-ins fall back to this

#ifndef BOOST_LIKELY
#define BOOST_LIKELY(x) (x)
#endif

#ifndef BOOST_UNLIKELY
#define BOOST_UNLIKELY(x) (x)
#endif

#include <cassert>
#include <cmath>

namespace salmon {

namespace math {
constexpr double LOG_0 = HUGE_VAL;
constexpr double LOG_1 = 0;
constexpr double LOG_ONEHALF = -0.69314718055994530941;
constexpr double LOG_ORPHAN_PROB = -2.30258509299404568401;
constexpr double EPSILON = 0.375e-10;
const double LOG_EPSILON = log(EPSILON);

inline double log(double v) { return (v > 0) ? std::log(v) : LOG_0; }
inline double exp(double v) { return std::exp(v); }
inline bool isLog0(double v) { return v == LOG_0; }
// Taken from https://github.com/adarob/eXpress/blob/master/src/main.h
inline bool approxEqual(double a, double b, double eps = EPSILON) {
  return std::abs(a - b) <= eps;
}

// Taken from https://github.com/adarob/eXpress/blob/master/src/main.h
inline double logAdd(double x, double y) {
  if (std::abs(x) == LOG_0) {
    return y;
  }
  if (std::abs(y) == LOG_0) {
    return x;
  }
  if (y > x) {
    std::swap(x, y);
  }
  double sum = x + std::log(1 + std::exp(y - x));
  return sum;
}

// Taken from https://github.com/adarob/eXpress/blob/master/src/main.h
inline double logSub(double x, double y) {
  if (std::abs(y) == LOG_0) {
    return x;
  }
  if (x <= y) {
    assert(std::fabs(x - y) < 1e-5);
    return LOG_0;
  }
  double diff = x + std::log(1 - std::exp(y - x));
  return diff;
}

} // namespace math

} // namespace salmon

#endif // SALMON_MATH_HPP
