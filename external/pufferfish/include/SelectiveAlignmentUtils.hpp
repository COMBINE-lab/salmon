#ifndef __SELECTIVE_ALIGNMENT_UTILS__
#define __SELECTIVE_ALIGNMENT_UTILS__

#include "PuffAligner.hpp"

namespace selective_alignment {
  namespace utils {


inline bool recoverOrphans(std::string& leftRead,
                    std::string& rightRead,
                    std::vector<pufferfish::util::MemCluster> &recoveredMemClusters,
                    std::vector<pufferfish::util::JointMems> &jointMemsList,
                    PuffAligner& puffaligner,
                    bool verbose) {
  using MateStatus = pufferfish::util::MateStatus;

  puffaligner.orphanRecoveryMemCollection.reserve(2 * jointMemsList.size() + 1);
  recoveredMemClusters.reserve(2 * jointMemsList.size() + 1);
  // if we recovered any orphans, then we discard all mappings
  // where we did not recover the mate.
  bool recoveredAny{false};
  for (auto& jointMem : jointMemsList) {
    auto& leftClust = jointMem.leftClust;
    auto& rightClust = jointMem.rightClust;
    auto tid = jointMem.tid;
    if (jointMem.mateStatus == MateStatus::PAIRED_END_LEFT) {
      bool recovered = puffaligner.recoverSingleOrphan(leftRead, rightRead, *leftClust, recoveredMemClusters, tid, true, verbose);
     if (recovered) {
       jointMem.rightClust = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
       jointMem.recovered = true;
       jointMem.mateStatus = MateStatus::PAIRED_END_PAIRED;
       recoveredAny = true;
    }
  } else if (jointMem.mateStatus == MateStatus::PAIRED_END_RIGHT) {
      bool recovered = puffaligner.recoverSingleOrphan(leftRead, rightRead, *rightClust, recoveredMemClusters, tid, false, verbose);
      if (recovered) {
        jointMem.leftClust = recoveredMemClusters.begin() + recoveredMemClusters.size() - 1;
        jointMem.recovered = true;
        jointMem.mateStatus = MateStatus::PAIRED_END_PAIRED;
        recoveredAny = true;
      }
    }
  }

  if (recoveredAny) {
    jointMemsList.erase(std::remove_if(jointMemsList.begin(), jointMemsList.end(),
                                      [](const pufferfish::util::JointMems& jm) { return !jm.recovered; }),
                       jointMemsList.end());
  }

  return recoveredAny;
}

  }
}

#endif // __SELECTIVE_ALIGNMENT_UTILS__
