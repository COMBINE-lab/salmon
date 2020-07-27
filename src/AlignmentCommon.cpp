#include "AlignmentCommon.hpp"

#include "ReadPair.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"
#include "SalmonStringUtils.hpp"

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

bool AlignmentCommon::computeErrorCount(bam_seq_t* read, Transcript& ref, ErrorCount& counts) {
  using namespace salmon::stringtools;

  auto         transcriptIdx = bam_pos(read);
  const size_t transcriptLen = ref.RefLength;
  int32_t readIdx{0};

  // if the read starts before the beginning of the transcript,
  // only consider the part overlapping the transcript
  if (transcriptIdx < 0) {
    readIdx = -transcriptIdx;
    transcriptIdx = 0;
  }
  // unsigned version of transcriptIdx
  size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

  uint32_t*      cigar       = bam_cigar(read);
  const uint32_t cigarLen    = bam_cigar_len(read);
  if (cigarLen == 0 || !cigar)
    return false;

  uint8_t*       qseq        = reinterpret_cast<uint8_t*>(bam_seq(read));
  const int32_t  readLen     = bam_seq_len(read);

  const strand readStrand = strand::forward;
  bool advanceInRead{false};
  bool advanceInReference{false};

  uint32_t errors = 0;
  uint32_t clips  = 0;

  for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {
    uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
    enum cigar_op op =
      static_cast<enum cigar_op>(cigar[cigarIdx] & BAM_CIGAR_MASK);

    size_t curReadBase = (BAM_CONSUME_SEQ(op)) ? samToTwoBit[bam_seqi(qseq, readIdx)] : 0;
    size_t curRefBase = (BAM_CONSUME_REF(op)) ? samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)] : 0;
    advanceInReference = false;

    for (size_t i = 0; i < opLen; ++i) {
      if (advanceInRead) {
        // Shouldn't happen!
        if (readIdx >= readLen) {
          if (logger_) {
            logger_->warn("CIGAR string for read [{}] "
                          "seems inconsistent. It refers to non-existant "
                          "positions in the read ({} >= {})!",
                          bam_name(read), readIdx, readLen);
            std::stringstream cigarStream;
            for (size_t j = 0; j < cigarLen; ++j) {
              uint32_t opLen = cigar[j] >> BAM_CIGAR_SHIFT;
              enum cigar_op op =
                static_cast<enum cigar_op>(cigar[j] & BAM_CIGAR_MASK);
              cigarStream << opLen << opToChr(op);
            }
            logger_->warn("CIGAR = {}", cigarStream.str());
          }
          return false;
        }

        curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        advanceInRead = false;
      }
      if (advanceInReference) {
        // Shouldn't happen!
        if (uTranscriptIdx >= transcriptLen) {
          if (logger_) {
            logger_->warn(
                          "CIGAR string for read [{}] "
                          "seems inconsistent. It refers to non-existant "
                          "positions in the reference! Transcript name "
                          "is {}, length is {}, id is {}. Read things refid is {}",
                          bam_name(read), ref.RefName, transcriptLen, ref.id,
                          bam_ref(read));
          }
          return false;
        }

        curRefBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
        advanceInReference = false;
      }
      setBasesFromCIGAROp_(op, curRefBase, curReadBase);
      if(curRefBase == ALN_SOFT_CLIP)
        ++clips;
      else if(curRefBase != curReadBase)
        ++errors;

      if (BAM_CONSUME_SEQ(op)) {
        ++readIdx;
        advanceInRead = true;
      }
      if (BAM_CONSUME_REF(op)) {
        ++uTranscriptIdx;
        advanceInReference = true;
      }
    }
  }

  counts.ims   = errors;
  counts.clips = clips;
  return true;
}
