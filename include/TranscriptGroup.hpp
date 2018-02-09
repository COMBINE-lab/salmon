#ifndef TRANSCRIPT_GROUP_HPP
#define TRANSCRIPT_GROUP_HPP

#include <boost/functional/hash.hpp>

#include <cstdint>
#include <vector>
#include "xxhash.h"

class TranscriptGroup {
public:
  TranscriptGroup();
  TranscriptGroup(std::vector<uint32_t> txpsIn);

  TranscriptGroup(std::vector<uint32_t> txpsIn, size_t hashIn);

  TranscriptGroup(TranscriptGroup&& other);
  TranscriptGroup(const TranscriptGroup& other);

  TranscriptGroup& operator=(const TranscriptGroup& other);
  TranscriptGroup& operator=(TranscriptGroup&& other);

  friend bool operator==(const TranscriptGroup& lhs,
                         const TranscriptGroup& rhs);

  void setValid(bool v) const;

  std::vector<uint32_t> txps;
  size_t hash;
  double totalMass;
  mutable bool valid;
};

bool operator==(const TranscriptGroup& lhs, const TranscriptGroup& rhs);

struct TranscriptGroupHasher {
  std::size_t operator()(const TranscriptGroup& k) const {
    return k.hash;
    /*
    std::size_t seed{0};
    for (auto e : k.txps) {
        boost::hash_combine(seed, e);
    }
    return seed;
    */
  }
};

#endif // TRANSCRIPT_GROUP_HPP
