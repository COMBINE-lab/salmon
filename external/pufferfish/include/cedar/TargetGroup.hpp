#ifndef TARGET_GROUP_HPP
#define TARGET_GROUP_HPP

#include <cstdint>
#include <cstddef>
#include <vector>
#include "xxhash.h"

class TargetGroup {
public:
  TargetGroup();
  TargetGroup(std::vector<uint32_t> tgtIn);

  TargetGroup(std::vector<uint32_t> tgtIn, size_t hashIn);

  TargetGroup(TargetGroup&& other);
  TargetGroup(const TargetGroup& other);

  TargetGroup& operator=(const TargetGroup& other);
  TargetGroup& operator=(TargetGroup&& other);

  friend bool operator==(const TargetGroup& lhs,
                         const TargetGroup& rhs);

  void setValid(bool v) const;

  std::vector<uint32_t> tgts;
  size_t hash;
  double totalMass;
  mutable bool valid;
};

bool operator==(const TargetGroup& lhs, const TargetGroup& rhs);

struct TargetGroupHasher{
  std::size_t operator()(const TargetGroup& k) const {
    return k.hash;
  }
};

#endif // TARGET GROUP 
