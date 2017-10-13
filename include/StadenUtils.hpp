#ifndef STADEN_UTILS
#define STADEN_UTILS

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

#include <cstdlib>

namespace staden {
namespace utils {
bam_seq_t* bam_init();
void bam_destroy(bam_seq_t* b);
} // namespace utils
} // namespace staden

#endif // STADEN_UTILS
