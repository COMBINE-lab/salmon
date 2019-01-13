#ifndef __SALMON_MAPPING_STATISTICS__
#define __SALMON_MAPPING_STATISTICS__

#include <atomic>

class MappingStatistics {
public:
  std::atomic<uint64_t> numOrphansRescued{0};
};

#endif // __SALMON_MAPPING_STATISTICS__
