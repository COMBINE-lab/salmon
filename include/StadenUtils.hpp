#ifndef STADEN_UTILS
#define STADEN_UTILS

extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
#undef max
#undef min
}

#include <cstdlib>

namespace staden {
    namespace utils {
        bam_seq_t* bam_init();
        void bam_destroy(bam_seq_t* b);
    }
}

#endif // STADEN_UTILS
