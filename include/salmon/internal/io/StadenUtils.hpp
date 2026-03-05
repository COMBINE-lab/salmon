#ifndef STADEN_UTILS
#define STADEN_UTILS

#include "salmon/internal/io/AlignmentIO.hpp"

namespace staden {
namespace utils {
bam_seq_t* bam_init();
void bam_destroy(bam_seq_t* b);
} // namespace utils
} // namespace staden

#endif // STADEN_UTILS
