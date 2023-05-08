#ifndef PUFFERFISH_TYPES_HPP
#define PUFFERFISH_TYPES_HPP

#include "CanonicalKmerIterator.hpp"

namespace pufferfish {
    namespace types {
        
using hasher_t = boomphf::SingleHashFunctor<uint64_t>;
using boophf_t = boomphf::mphf<uint64_t, hasher_t>;
using EqClassID = uint32_t;
using EqClassLabel = std::vector<uint32_t>;
using CanonicalKmerIterator = pufferfish::CanonicalKmerIterator ;
 
    }
}
#endif //PUFFERFISH_TYPES_HPP