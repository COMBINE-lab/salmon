#ifndef __SALMON_MAPPING_STATISTICS__
#define __SALMON_MAPPING_STATISTICS__

#include <atomic>

// class to collect statistics we may want about mapping / alignment
class MappingStatistics {
public:
  std::atomic<uint64_t> numOrphansRescued{0};
  // The number of mappings that were filtered due to alignment score
  std::atomic<uint64_t> numMappingsFiltered{0};
  std::atomic<uint64_t> numFragmentsFiltered{0};
  std::atomic<uint64_t> numDecoyFragments{0};
  std::atomic<uint64_t> numDovetails{0};
};

#endif // __SALMON_MAPPING_STATISTICS__
