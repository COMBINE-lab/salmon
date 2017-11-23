#ifndef FASTX_PARSER_THREAD_UTILS_HPP
#define FASTX_PARSER_THREAD_UTILS_HPP

#include <cassert>
#include <chrono>
#include <pthread.h>
#include <random>
#include <thread>

// Most of this code is taken directly from
// https://github.com/geidav/spinlocks-bench/blob/master/os.hpp. However, things
// may be renamed, modified, or randomly mangled over time.
#define ALWAYS_INLINE inline __attribute__((__always_inline__))

namespace fastx_parser {
namespace thread_utils {

static const constexpr size_t MIN_BACKOFF_ITERS = 32;
static const size_t MAX_BACKOFF_ITERS = 1024;

ALWAYS_INLINE static void cpuRelax() { asm("pause"); }

ALWAYS_INLINE void yieldSleep() {
  using namespace std::chrono;
  std::chrono::microseconds ytime(500);
  std::this_thread::sleep_for(ytime);
}

ALWAYS_INLINE void backoffExp(size_t& curMaxIters) {
  thread_local std::uniform_int_distribution<size_t> dist;
  thread_local std::minstd_rand gen(std::random_device{}());
  const size_t spinIters =
      dist(gen, decltype(dist)::param_type{0, curMaxIters});
  curMaxIters = std::min(2 * curMaxIters, MAX_BACKOFF_ITERS);
  for (size_t i = 0; i < spinIters; i++) {
    cpuRelax();
  }
}

ALWAYS_INLINE void backoffOrYield(size_t& curMaxDelay) {
  if (curMaxDelay >= MAX_BACKOFF_ITERS) {
    yieldSleep();
    curMaxDelay = MIN_BACKOFF_ITERS;
  }
  backoffExp(curMaxDelay);
}

} // namespace thread_utils
} // namespace fastx_parser

#endif // FASTX_PARSER_THREAD_UTILS_HPP
