#ifndef UTILITY_FUNCTIONS_HPP
#define UTILITY_FUNCTIONS_HPP

#include "SalmonUtils.hpp"
#include <limits>

// from
// http://stackoverflow.com/questions/17719674/c11-fast-constexpr-integer-powers
constexpr int64_t constExprPow(int64_t base, unsigned int exp,
                               int64_t result = 1) {
  return exp < 1 ? result
                 : constExprPow(base * base, exp / 2,
                                (exp % 2) ? result * base : result);
}

inline std::string kmerForIndex(uint32_t idx, uint32_t K) {
  std::string kmer(K, 'X');
  // The number of bits we need to shift the
  // current mask to the left.
  uint32_t pos{0};
  for (int32_t i = K - 1; i >= 0; --i) {
    uint8_t c = (idx >> pos) & 0x3;
    switch (c) {
    case 0:
      kmer[i] = 'A';
      break;
    case 1:
      kmer[i] = 'C';
      break;
    case 2:
      kmer[i] = 'G';
      break;
    case 3:
      kmer[i] = 'T';
      break;
    default:
      break;
    }
    pos += 2;
  }
  return kmer;
}

inline uint32_t nextKmerIndex(uint32_t idx, char n, uint32_t K,
                              salmon::utils::Direction dir) {
  using salmon::utils::Direction;
  if (dir == Direction::REVERSE or dir == Direction::REVERSE_COMPLEMENT) {
    // drop the leftmost character, and replace it with the complement of the
    // new one.
    idx = idx >> 2;
    switch (n) {
    case 'A':
    case 'a':
      //  n='T';
      // complement is 'T';
      idx = idx | (3 << 2 * (K - 1));
      break;
    case 'C':
    case 'c':
      //  n='G';
      // complement is 'G';
      idx = idx | (2 << 2 * (K - 1));
      break;
    case 'g':
    case 'G':
      // n='C';
      // complement is 'C';
      idx = idx | (1 << 2 * (K - 1));
      break;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      // n='A';
      // complement is 'A';
      break;
    }
    return idx;
  } else {
    // drop the rightmost character and replace it with the new one.
    idx = idx << 2;
    switch (n) {
    case 'A':
    case 'a':
      break;
    case 'C':
    case 'c':
      idx = idx + 1;
      break;
    case 'G':
    case 'g':
      idx = idx + 2;
      break;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      idx = idx + 3;
      break;
    }
    // Clear the top 32 - 2*K bits.
    uint32_t clearShift = (32 - 2 * K);
    return idx & (0xFFFFFFFF >> clearShift);
  }
}

inline uint32_t indexForKmer(const char* s, uint32_t K,
                             salmon::utils::Direction dir) {
  using salmon::utils::Direction;
  // The index we'll return
  uint32_t idx{0};
  int32_t SK = static_cast<int32_t>(K);
  // The number of bits we need to shift the
  // current mask to the left.
  if (dir == Direction::FORWARD) {
    for (int32_t i = 0; i < SK; ++i) {
      switch (s[i]) {
      case 'A':
      case 'a':
        break;
      case 'C':
      case 'c':
        idx += 1;
        break;
      case 'G':
      case 'g':
        idx += 2;
        break;
      case 'T':
      case 't':
      case 'U':
      case 'u':
        idx += 3;
        break;
      default:
        return std::numeric_limits<uint32_t>::max();
      }
      if (i < SK - 1) {
        idx = idx << 2;
      }
    }
  } else {
    for (int32_t i = SK - 1; i >= 0; i--) {
      switch (s[i]) {
      case 'T':
      case 't':
      case 'u':
      case 'U':
        break;
      case 'C':
      case 'c':
        idx += 2;
        break;
      case 'G':
      case 'g':
        idx += 1;
        break;
      case 'A':
      case 'a':
        idx += 3;
        break;
      default:
        return std::numeric_limits<uint32_t>::max();
      }
      if (i > 0) {
        idx = idx << 2;
      }
    }
  }
  return idx;
}

// Range iterator to fill vector upon construction (E.g., std::vector<int> v(range{0}, range{5}) -> [0, 1, 2, 3, 4]).
struct range : std::iterator<std::forward_iterator_tag,int> {
  int x;
  range(int x_) : x(x_) { }
  range& operator++() { ++x; return *this; }
  int operator*() const { return x; }
  int operator-(const range& rhs) const { return x - rhs.x; }
  bool operator!=(const range& rhs) const { return rhs.x != x; }
};

// Permute in place the elements of ary according to the permutation
// vector. permutation is destroyed in the process (all elements set
// to -1). The elements of permutation must be signed integers.
template<typename V1, typename V2, typename V3>
void permute(V1& permutation, V2& ary1, V3& ary2) {
  for(size_t i = 0; i < permutation.size(); ++i) {
    if(permutation[i] < 0) continue;
    int cur = i;
    int next = permutation[cur];
    permutation[cur] = -1;
    while(next != i) {
      std::swap(ary1[cur], ary1[next]);
      std::swap(ary2[cur], ary2[next]);
      cur = next;
      next = permutation[cur];
      permutation[cur] = -1;
    }
  }
}

#endif // UTILITY_FUNCTIONS_HPP
