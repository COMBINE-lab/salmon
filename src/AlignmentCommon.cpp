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

void AlignmentCommon::setBasesFromCIGAROp_(enum cigar_op op,
                                           size_t& curRefBase,
                                           size_t& curReadBase) {
  switch (op) {
  case BAM_UNKNOWN:
    std::cerr << "ENCOUNTERED UNKNOWN SYMBOL IN CIGAR STRING!" << std::endl;
    break;
  case BAM_CMATCH:
    // do nothing
    break;
  case BAM_CBASE_MATCH:
    // do nothing
    break;
  case BAM_CBASE_MISMATCH:
    // do nothing
    break;
  case BAM_CINS:
    curRefBase = ALN_DASH;
    break;
  case BAM_CDEL:
    curReadBase = ALN_DASH;
    break;
  case BAM_CREF_SKIP:
    curReadBase = ALN_REF_SKIP;
    break;
  case BAM_CSOFT_CLIP:
    curRefBase = ALN_SOFT_CLIP;
    break;
  case BAM_CHARD_CLIP:
    curRefBase = ALN_HARD_CLIP;
    curReadBase = ALN_HARD_CLIP;
    break;
  case BAM_CPAD:
    curRefBase = ALN_PAD;
    curReadBase = ALN_PAD;
    break;
  default:
    std::cerr << "ENCOUNTERED UNKNOWN (non -1) CIGAR OP : (" << op << ")!" << std::endl;
    break;
  }
}

char AlignmentCommon::opToChr(enum cigar_op op) {
  switch (op) {
  case BAM_UNKNOWN:
    std::cerr << "ENCOUNTERED UNKNOWN SYMBOL IN CIGAR STRING!\n";
    break;
  case BAM_CMATCH:
    // do nothing
    return 'M';
    break;
  case BAM_CBASE_MATCH:
    // do nothing
    return 'M';
    break;
  case BAM_CBASE_MISMATCH:
    // do nothing
    return 'M';
    break;
  case BAM_CINS:
    return 'I';
    break;
  case BAM_CDEL:
    return 'D';
    break;
  case BAM_CREF_SKIP:
    return 'S';
    break;
  case BAM_CSOFT_CLIP:
    return 'c';
    break;
  case BAM_CHARD_CLIP:
    return 'C';
    break;
  case BAM_CPAD:
    return 'P';
    break;
  }
  return 'X';
}
