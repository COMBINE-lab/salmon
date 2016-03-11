#ifndef READ_KMER_DIST_HPP
#define READ_KMER_DIST_HPP

#include <fstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <cstdint>
#include "UtilityFunctions.hpp"
#include "SalmonUtils.hpp"

template <uint32_t K, typename CountT = uint32_t>
class ReadKmerDist {
  public:
    std::array<CountT, constExprPow(4,K)> counts;

    ReadKmerDist() {
      // set a pseudo-count of 1
      for (size_t i = 0; i < counts.size(); ++i) {
	counts[i] = 1;
      }
    }

    inline constexpr uint32_t getK() { return K; }

    inline uint64_t totalCount() {
      CountT c{0};
      for (auto const& rc : counts) { c += rc; }
      return c;
    }

    // update the k-mer context for the hit at position p.
    // The underlying transcript is from [start, end)
    inline bool update(const char* start, const char *p, const char *end,
	salmon::utils::Direction dir) {
      using salmon::utils::Direction;
      int posBeforeHit = 2;
      int posAfterHit = 3;
      bool success{false};
      switch (dir) {
	case Direction::FORWARD :
	  {
	    // If we can fit the window before and after the read
	    if ((p - start) >= posBeforeHit and
		((p - posBeforeHit + K) < end) ) {
	      p -= posBeforeHit;
	      // If the read matches in the forward direction, we take
	      // the RC sequence.
	      auto idx = indexForKmer(p, K, Direction::REVERSE_COMPLEMENT);
	      if (idx > counts.size()) { return false; }
	      counts[idx]++;
	      success = true;
	    }
	  }
	  break;
	case Direction::REVERSE_COMPLEMENT :
	  {
	    if ((p - start) >= posAfterHit and
		((p - posAfterHit + K) < end) ) {
	      p -= posAfterHit;
	      auto idx = indexForKmer(p, K, Direction::FORWARD);
	      if (idx > counts.size()) { return false; }
	      counts[idx]++;
	      success = true;
	    }
	  }
	  break;
	default:
	  break;
      }
      return success;
    }

};
#endif // READ_KMER_DIST_HPP 
