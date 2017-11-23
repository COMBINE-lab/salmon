#ifndef MICROOPTS_HPP
#define MICROOPTS_HPP

/** Based on Likely.h from FB's Folly library **/

#undef LIKELY
#undef UNLIKELY

#if defined(__GUNC__) && __GNUC__ >= 4
#define LIKELY(X) (__builtin_expect((x), 1))
#define UNLIKELY(X) (__builtin_expect((x), 0))
#elif defined(__clang__)

#endif // MICROOPTS_HPP
