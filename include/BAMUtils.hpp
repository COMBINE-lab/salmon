#ifndef __SALMON_BAM_UTILS__
#define __SALMON_BAM_UTILS__

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

#include <unordered_map>
#include <string>

namespace salmon {
  namespace bam_utils {

    enum class AlignerDetails : uint8_t
      {
       BOWTIE2=1,     // a standard bowtie 2 configuration
       BOWTIE2_LOCAL, // using --local flag
       STAR,          // STAR aligner
       BWA_MEM,       // BWA-MEM aligner
       MINIMAP2,      // minimap2 aligner
       RAPMAP,        // rapmap aligner/mapper
       PUFFERFISH,    // pufferfish / puffalign aligner
       UNKNOWN
      };

    std::string to_string(AlignerDetails a);
    AlignerDetails inferAlignerFromHeader(SAM_hdr* header);
  }
}

#endif // __SALMON_BAM_UTILS__
