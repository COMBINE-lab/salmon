#include <cstdint>
#include <vector>

#include "cedar/TargetGroup.hpp"
#include "xxhash.h"

TargetGroup::TargetGroup() : hash(0) {}

TargetGroup::TargetGroup(std::vector<uint32_t> tgtsIn)
    : tgts(tgtsIn), valid(true) {
  size_t seed{0};
  hash = XXH64(static_cast<void*>(tgts.data()), tgts.size() * sizeof(uint32_t),
               seed);
}

TargetGroup::TargetGroup(std::vector<uint32_t> tgtsIn, size_t hashIn)
    : tgts(tgtsIn), hash(hashIn), valid(true) {}

TargetGroup::TargetGroup(const TargetGroup& other) {
  tgts = other.tgts;
  hash = other.hash;
  valid = other.valid;
}

TargetGroup& TargetGroup::operator=(const TargetGroup& other) {
  tgts = other.tgts;
  hash = other.hash;
  valid = other.valid;
  return *this;
}

TargetGroup::TargetGroup(TargetGroup&& other) {
  tgts = std::move(other.tgts);
  hash = other.hash;
  valid = other.valid;
}

void TargetGroup::setValid(bool b) const { valid = b; }

TargetGroup& TargetGroup::operator=(TargetGroup&& other) {
  tgts = std::move(other.tgts);
  hash = other.hash;
  valid = other.valid;
  return *this;
}

bool operator==(const TargetGroup& lhs, const TargetGroup& rhs) {
  return lhs.tgts == rhs.tgts;
};
