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

bool AlignmentCommon::computeErrorCount(bam_seq_t* read, bam_seq_t* primary, Transcript& ref, ErrorCount& counts,
                                        const char* src) {
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

  const bool    usePrimary = bam_seq_len(read) == 0 && primary != nullptr;
  uint8_t*      qseq       = reinterpret_cast<uint8_t*>(!usePrimary ? bam_seq(read) : bam_seq(primary));
  const int32_t readLen    = !usePrimary ? bam_seq_len(read) : bam_seq_len(primary);
  // If using the primary sequence and the reverse flag is different
  // between the current sequence and the primary, then we use the
  // reverse strand.
  const strand readStrand = usePrimary && (bam_flag(primary) & BAM_FREVERSE) != (bam_flag(read) & BAM_FREVERSE)
    ? strand::reverse : strand::forward;
  counts.clear();

  // Go through every operation in the CIGAR string as long as we
  // still point within the sequences. Note that here we allow readIdx
  // == readLen. This happens when there is a hard clip (H) as the
  // last operation and it is not an error. Similarly, it is possible
  // to have uTranscriptIdx == transcriptLen (read that extend beyond
  // the end of the reference) and it is not an error if the following
  // operations do not consuming reference sequence.
  for (uint32_t cigarIdx = 0; cigarIdx < cigarLen && readIdx <= readLen && uTranscriptIdx <= transcriptLen; ++cigarIdx) {
    uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
    enum cigar_op op =
      static_cast<enum cigar_op>(cigar[cigarIdx] & BAM_CIGAR_MASK);

    // Only "matching" (meaning match or mismatch need to look at each
    // base of query and reference sequences. Otherwise, advance block by block
    switch(op) {
    case BAM_CSOFT_CLIP:
      readIdx          += opLen;
      counts.sclips_   += opLen;
      if(cigarIdx == 0)
        counts.fclips_ += opLen;
      else
        counts.bclips_ += opLen;
      continue;

    case BAM_CDEL:
      counts.deletions_ += opLen; // Fallthrough
    case BAM_CREF_SKIP:
      uTranscriptIdx    += opLen;
      continue;

    case BAM_CBASE_MISMATCH:
      counts.mismatches_ += opLen;
      readIdx            += opLen;
      uTranscriptIdx     += opLen;
      continue;

    case BAM_CBASE_MATCH:
      counts.matches_ += opLen;
      readIdx         += opLen;
      uTranscriptIdx  += opLen;
      continue;

    case BAM_CINS:
      counts.insertions_ += opLen;
      readIdx            += opLen;
      continue;

    case BAM_CHARD_CLIP:
      counts.hclips_   += opLen;
      if(cigarIdx == 0)
        counts.fclips_ += opLen;
      else
        counts.bclips_ += opLen;      // Fallthrough
    case BAM_CPAD:
      continue;

    case BAM_CMATCH:
      break; // Examine base by base

    case BAM_UNKNOWN:
      if(logger_)
        logger_->warn("in {} found unknown CIGAR operation for read {}",
                      src, bam_name(read));
      return false;

    default:
      if(logger_)
        logger_-> warn("in {} found a bug, parsing read {}", src, bam_name(read));
      return false;
    }

    // "Match". Count mismatches
    for (size_t i = 0; i < opLen && readIdx < readLen && uTranscriptIdx < transcriptLen; ++i) {
      const auto curReadBase  =
        readStrand           == strand::forward
        ? samToTwoBit[bam_seqi(qseq, readIdx)]
        : 3 - samToTwoBit[bam_seqi(qseq, readLen - readIdx - 1)];
      const auto curRefBase   = samToTwoBit[ref.baseAt(uTranscriptIdx)];
      const int  isMatch      = curReadBase == curRefBase;
      counts.matches_        += isMatch;
      counts.mismatches_     += !isMatch;

      ++readIdx;
      ++uTranscriptIdx;
    }
  }

  if (readIdx > readLen) {
    if (logger_) {
      std::stringstream cigarStream;
      for (size_t j = 0; j < cigarLen; ++j) {
        uint32_t opLen = cigar[j] >> BAM_CIGAR_SHIFT;
        enum cigar_op op =
          static_cast<enum cigar_op>(cigar[j] & BAM_CIGAR_MASK);
        cigarStream << opLen << opToChr(op);
      }
      logger_->warn("in {} CIGAR [{}] string for read [{}] "
                    "seems inconsistent. It refers to non-existant "
                    "positions in the read ({} >= {})! {}",
                    src, cigarStream.str(), bam_name(read), readIdx, readLen, usePrimary);
    }
    return false;
  }
  if (uTranscriptIdx > transcriptLen) {
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
  return true;
}
