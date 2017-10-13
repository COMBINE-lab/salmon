#include "StadenUtils.hpp"

namespace staden {
namespace utils {

bam_seq_t* bam_init() {
  return reinterpret_cast<bam_seq_t*>(calloc(1, sizeof(bam_seq_t)));
}

void bam_destroy(bam_seq_t* b) {
  if (b == 0) {
    return;
  }
  free(b);
}

} // namespace utils
} // namespace staden
