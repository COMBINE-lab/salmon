#ifndef __SALMON_MAPPING_STATISTICS__
#define __SALMON_MAPPING_STATISTICS__

#include <atomic>

// class to collect statistics we may want about mapping / alignment
class MappingStatistics {
public:
  std::atomic<uint64_t> numOrphansRescued{0};
};

#endif // __SALMON_MAPPING_STATISTICS__
