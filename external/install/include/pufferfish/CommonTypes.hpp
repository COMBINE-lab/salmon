#ifndef __PUFFERFISH_COMMON_TYPES_HPP__
#define __PUFFERFISH_COMMON_TYPES_HPP__

#include "core/range.hpp"
#include "itlib/small_vector.hpp"
#include "parallel_hashmap/phmap.h"
#include "Util.hpp"

namespace pufferfish {
  namespace common_types {
    using ReferenceID = size_t;
    template <typename T>
    using IterRange = core::range<typename T::iterator>;
    using RefMemMapT = pufferfish::util::CachedVectorMap<std::pair<pufferfish::common_types::ReferenceID, bool>, pufferfish::util::MemInfo, pufferfish::util::pair_hash>;
  }
}

#endif // __PUFFERFISH_COMMON_TYPES_HPP__
