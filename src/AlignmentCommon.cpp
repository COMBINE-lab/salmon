#include "AlignmentCommon.hpp"

#include "ReadPair.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"

bool AlignmentCommon::hasIndel(ReadPair& hit) {
  if (!hit.isPaired())
    return hasIndel(hit.read1);
  return hasIndel(hit.read1) or hasIndel(hit.read2);
}

bool AlignmentCommon::hasIndel(UnpairedRead& hit) {
  return hasIndel(hit.read);
}


bool AlignmentCommon::hasIndel(bam_seq_t* read)  {
  uint32_t* cigar = bam_cigar(read);
  uint32_t cigarLen = bam_cigar_len(read);

  for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {
    uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
    enum cigar_op op =
        static_cast<enum cigar_op>(cigar[cigarIdx] & BAM_CIGAR_MASK);

    switch (op) {
    case BAM_CINS:
      return true;
    case BAM_CDEL:
      return true;
    default:
      break;
    }
  }
  return false;
}
