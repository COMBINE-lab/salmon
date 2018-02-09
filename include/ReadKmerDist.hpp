#ifndef READ_KMER_DIST_HPP
#define READ_KMER_DIST_HPP

#include "SalmonUtils.hpp"
#include "UtilityFunctions.hpp"
#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <vector>

template <uint32_t K, typename CountT = uint32_t> class ReadKmerDist {
public:
  std::array<CountT, constExprPow(4, K)> counts;
  std::unordered_map<std::string, CountT> strCounts;
  std::map<char, uint32_t> startMap;

  ReadKmerDist(bool pseudoCount = true) {
    // set a pseudo-count of 1
    if (pseudoCount) {
      for (size_t i = 0; i < counts.size(); ++i) {
        counts[i] = 1;
      }
    }
  }

  inline constexpr uint32_t getK() const { return K; }

  inline uint64_t totalCount() {
    CountT c{0};
    for (auto const& rc : counts) {
      c += rc;
    }
    return c;
  }

  // update the k-mer context for the hit at position p.
  // The underlying transcript is from [start, end)
  inline bool update(const char* start, const char* p, const char* end,
                     salmon::utils::Direction dir) {
    using salmon::utils::Direction;
    // int posBeforeHit = 3;
    // int posAfterHit = 2;
    int posBeforeHit = 4;
    int posAfterHit = 3;
    bool success{false};
    bool contextExists{false};

    if (dir == Direction::FORWARD) {
      // If we can fit the window before and after the read
      if ((p - start) >= posBeforeHit and ((p - posBeforeHit + K) < end)) {
        p -= posBeforeHit;
        contextExists = true;
      }
    } else if (dir == Direction::REVERSE_COMPLEMENT) {
      if ((p - start) >= posAfterHit and ((p - posAfterHit + K) < end)) {
        p -= posAfterHit;
        contextExists = true;
      }
    }

    auto idx = contextExists ? indexForKmer(p, K, dir)
                             : std::numeric_limits<uint32_t>::max();

    if (idx < counts.size()) {
      counts[idx]++;
      success = true;
    }
    return success;
  }
};
#endif // READ_KMER_DIST_HPP
