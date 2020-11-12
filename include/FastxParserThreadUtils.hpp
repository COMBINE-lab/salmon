#ifndef FASTX_PARSER_THREAD_UTILS_HPP
#define FASTX_PARSER_THREAD_UTILS_HPP

#include <cassert>
#include <chrono>
#include <pthread.h>
#include <random>
#include <thread>
#if defined(__SSE2__)
#include "simde/x86/sse2.h"
//#include <xmmintrin.h> // _mm_pause
#endif

// Most of this code is taken directly from
// https://github.com/geidav/spinlocks-bench/blob/master/os.hpp. However, things
// may be renamed, modified, or randomly mangled over time.
#define ALWAYS_INLINE inline __attribute__((__always_inline__))

namespace fastx_parser {
namespace thread_utils {

static const constexpr size_t MIN_BACKOFF_ITERS = 32;
static const size_t MAX_BACKOFF_ITERS = 1024;

ALWAYS_INLINE static void cpuRelax() {
#if defined(__SSE2__)  // AMD and Intel
  simde_mm_pause();
#elif defined(__i386__) || defined(__x86_64__)
  asm volatile("pause");
#elif defined(__aarch64__)
  asm volatile("wfe");
#elif defined(__armel__) || defined(__ARMEL__)
  asm volatile ("nop" ::: "memory");
#elif defined(__arm__) || defined(__aarch64__)
  __asm__ __volatile__ ("yield" ::: "memory");
#elif defined(__ia64__)  // IA64
  __asm__ __volatile__ ("hint @pause");
#elif defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
   __asm__ __volatile__ ("or 27,27,27" ::: "memory");
#else  // everything else.
   asm volatile ("nop" ::: "memory");
#endif
}

ALWAYS_INLINE void yieldSleep() {
  using namespace std::chrono;
  std::chrono::microseconds ytime(500);
  std::this_thread::sleep_for(ytime);
}

ALWAYS_INLINE void backoffExp(size_t& curMaxIters) {
  thread_local std::uniform_int_distribution<size_t> dist;

  // see : https://github.com/coryan/google-cloud-cpp-common/blob/a6e7b6b362d72451d6dc1fec5bc7643693dbea96/google/cloud/internal/random.cc
  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    thread_local std::random_device rd("/dev/urandom");
  #else
    thread_local std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128

  thread_local std::minstd_rand gen(rd());
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
