#include <cstdint>
#include <vector>

#include "SalmonMath.hpp"
#include "TranscriptGroup.hpp"
#include "xxhash.h"

TranscriptGroup::TranscriptGroup() : hash(0) {}

TranscriptGroup::TranscriptGroup(std::vector<uint32_t> txpsIn)
    : txps(txpsIn), valid(true) {
  size_t seed{0};
  hash = XXH64(static_cast<void*>(txps.data()), txps.size() * sizeof(uint32_t),
               seed);
}

TranscriptGroup::TranscriptGroup(std::vector<uint32_t> txpsIn, size_t hashIn)
    : txps(txpsIn), hash(hashIn), valid(true) {}

TranscriptGroup::TranscriptGroup(const TranscriptGroup& other) {
  txps = other.txps;
  hash = other.hash;
  valid = other.valid;
}

TranscriptGroup& TranscriptGroup::operator=(const TranscriptGroup& other) {
  txps = other.txps;
  hash = other.hash;
  valid = other.valid;
  return *this;
}

TranscriptGroup::TranscriptGroup(TranscriptGroup&& other) {
  txps = std::move(other.txps);
  hash = other.hash;
  valid = other.valid;
}

void TranscriptGroup::setValid(bool b) const { valid = b; }

TranscriptGroup& TranscriptGroup::operator=(TranscriptGroup&& other) {
  txps = std::move(other.txps);
  hash = other.hash;
  valid = other.valid;
  return *this;
}

bool operator==(const TranscriptGroup& lhs, const TranscriptGroup& rhs) {
  return lhs.txps == rhs.txps;
};
