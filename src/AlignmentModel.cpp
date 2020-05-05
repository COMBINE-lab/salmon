#include <cstdio>
#include <iostream>
#include <map>
#include <tuple>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY

#include "AlignmentModel.hpp"
#include "ReadPair.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"

AlignmentModel::AlignmentModel(double alpha, uint32_t readBins)
    : transitionProbsLeft_(readBins), transitionProbsRight_(readBins),
      isEnabled_(true), readBins_(readBins), burnedIn_(false) {

  for (size_t i = 0; i < readBins; ++i) {
    transitionProbsLeft_[i] = std::move(AtomicMatrix<double>(
        numAlignmentStates(), numAlignmentStates(), alpha));
    transitionProbsRight_[i] = std::move(AtomicMatrix<double>(
        numAlignmentStates(), numAlignmentStates(), alpha));
  }
}

bool AlignmentModel::burnedIn() { return burnedIn_; }
void AlignmentModel::burnedIn(bool burnedIn) { burnedIn_ = burnedIn; }

void AlignmentModel::setLogger(std::shared_ptr<spdlog::logger> logger) {
  logger_ = logger;
}

bool AlignmentModel::hasLogger() { return (logger_) ? true : false; }

bool AlignmentModel::hasIndel(ReadPair& hit) {
  if (!hit.isPaired()) {
    return hasIndel(hit.read1);
  }
  return hasIndel(hit.read1) or hasIndel(hit.read2);
}

bool AlignmentModel::hasIndel(UnpairedRead& hit) { return hasIndel(hit.read); }

bool AlignmentModel::hasIndel(bam_seq_t* read) {
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

/*
inline void AlignmentModel::setBasesFromCIGAROp_(enum cigar_op op, size_t&
curRefBase, size_t& curReadBase, std::stringstream& readStr, std::stringstream&
matchStr, std::stringstream& refStr) { using salmon::stringtools::twoBitToChar;
    switch (op) {
        case BAM_UNKNOWN:
            std::cerr << "ENCOUNTERED UNKNOWN SYMBOL IN CIGAR STRING!\n";
            break;
        case BAM_CMATCH:
            readStr << twoBitToChar[curReadBase];
            matchStr << ((curReadBase == curRefBase) ? ' ' : 'X');
            refStr << twoBitToChar[curRefBase];
            // do nothing
            break;
        case BAM_CBASE_MATCH:
            readStr << twoBitToChar[curReadBase];
            matchStr << ' ';
            refStr << twoBitToChar[curRefBase];
            // do nothing
            break;
        case BAM_CBASE_MISMATCH:
            readStr << twoBitToChar[curReadBase];
            matchStr << 'X';
            refStr << twoBitToChar[curRefBase];
            // do nothing
            break;
        case BAM_CINS:
            readStr << twoBitToChar[curReadBase];
            matchStr << ' ';
            refStr << '-';
            curRefBase = ALN_DASH;
            break;
        case BAM_CDEL:
            readStr << '-';
            matchStr << ' ';
            refStr << twoBitToChar[curRefBase];
            curReadBase = ALN_DASH;
            break;
        case BAM_CREF_SKIP:
            readStr << 'N';
            matchStr << ' ';
            refStr << twoBitToChar[curRefBase];
            curReadBase = ALN_REF_SKIP;
            break;
        case BAM_CSOFT_CLIP:
            readStr << twoBitToChar[curReadBase];
            matchStr << ' ';
            refStr << 'S';
            curRefBase = ALN_SOFT_CLIP;
            break;
        case BAM_CHARD_CLIP:
            readStr << 'H';
            matchStr << ' ';
            refStr << 'H';
            curRefBase = ALN_HARD_CLIP;
            curReadBase = ALN_HARD_CLIP;
            break;
        case BAM_CPAD:
            readStr << 'P';
            matchStr << ' ';
            refStr << 'P';
            curRefBase = ALN_PAD;
            curReadBase = ALN_PAD;
            break;
    }
}
*/

inline void AlignmentModel::setBasesFromCIGAROp_(enum cigar_op op,
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

char opToChr(enum cigar_op op) {
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

AlignmentModel::AlnModelProb AlignmentModel::logLikelihood(
    bam_seq_t* read, Transcript& ref,
    std::vector<AtomicMatrix<double>>& transitionProbs) {
  using namespace salmon::stringtools;
  bool useQual{false};
  size_t readIdx{0};
  auto transcriptIdx = bam_pos(read);
  size_t transcriptLen = ref.RefLength;
  // if the read starts before the beginning of the transcript,
  // only consider the part overlapping the transcript
  if (transcriptIdx < 0) {
    readIdx = -transcriptIdx;
    transcriptIdx = 0;
  }

  // unsigned version of transcriptIdx
  size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

  if (uTranscriptIdx >= transcriptLen) {
    std::lock_guard<std::mutex> l(outputMutex_);
    std::cerr << "transcript index = " << uTranscriptIdx
              << ", transcript length = " << transcriptLen << "\n";
    return {salmon::math::LOG_0, salmon::math::LOG_1};
  }

  // std::stringstream readStream, matchStream, refStream;

  uint32_t* cigar = bam_cigar(read);
  uint32_t cigarLen = bam_cigar_len(read);
  uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));
  uint8_t* qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));
  int32_t readLen = bam_seq_len(read);

  if (cigarLen == 0 or !cigar) {
    return {salmon::math::LOG_EPSILON, salmon::math::LOG_1};
  }

  salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
  // foreground likelihood
  double logLike = salmon::math::LOG_1;
  // background likelihood
  double bgLogLike = salmon::math::LOG_1;

  bool advanceInRead{false};
  bool advanceInReference{false};
  uint32_t readPosBin{0};
  uint32_t cigarIdx{0};
  uint32_t prevStateIdx{startStateIdx};
  uint32_t curStateIdx{0};
  double invLen = static_cast<double>(readBins_) / bam_seq_len(read);

  for (uint32_t cigarIdx = 0; cigarIdx < cigarLen; ++cigarIdx) {
    uint32_t opLen = cigar[cigarIdx] >> BAM_CIGAR_SHIFT;
    enum cigar_op op =
        static_cast<enum cigar_op>(cigar[cigarIdx] & BAM_CIGAR_MASK);
    size_t curReadBase = (BAM_CONSUME_SEQ(op)) ? samToTwoBit[bam_seqi(qseq, readIdx)] : 0;
    size_t curRefBase = (BAM_CONSUME_REF(op)) ? samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)] : 0;
    advanceInRead = false;
    advanceInReference = false;

    for (size_t i = 0; i < opLen; ++i) {
      if (advanceInRead) {
        // Shouldn't happen!
        if (readIdx >= static_cast<decltype(readIdx)>(readLen)) {
          if (logger_) {
            logger_->warn("(in logLikelihood()) CIGAR string for read [{}] "
                          "seems inconsistent. It refers to non-existant "
                          "positions in the read!",
                          bam_name(read));
            std::stringstream cigarStream;
            for (size_t j = 0; j < cigarLen; ++j) {
              uint32_t opLen = cigar[j] >> BAM_CIGAR_SHIFT;
              enum cigar_op op =
                  static_cast<enum cigar_op>(cigar[j] & BAM_CIGAR_MASK);
              cigarStream << opLen << opToChr(op);
            }
            logger_->warn("(in logLikelihood()) CIGAR = {}", cigarStream.str());
          }
          return {logLike, bgLogLike};
        }
        curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
        readPosBin = static_cast<uint32_t>((readIdx * invLen));
        advanceInRead = false;
      }
      if (advanceInReference) {
        // Shouldn't happen!
        if (uTranscriptIdx >= transcriptLen) {
          if (logger_) {
            logger_->warn(
                "(in logLikelihood()) CIGAR string for read [{}] "
                "seems inconsistent. It refers to non-existant "
                "positions in the reference! Transcript name "
                "is {}, length is {}, id is {}. Read things refid is {}",
                bam_name(read), ref.RefName, transcriptLen, ref.id,
                bam_ref(read));
          }
          return {logLike, bgLogLike};
        }
        curRefBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
        advanceInReference = false;
      }

      setBasesFromCIGAROp_(
          op, curRefBase, curReadBase); //, readStream, matchStream, refStream);

      curStateIdx = curRefBase * numStates + curReadBase;
      double tp = transitionProbs[readPosBin](prevStateIdx, curStateIdx);
      logLike += tp;
      bgLogLike += transitionProbs[readPosBin](0, 0);
      prevStateIdx = curStateIdx;
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

  /*
  {
      std::lock_guard<std::mutex> l(outputMutex_);
      std::cerr << "\n\nread:   " << readStream.str() << "\n";
      std::cerr << "        " << matchStream.str() << "\n";
      std::cerr << "ref:    " << refStream.str() << "\n";
      //std::cerr << "states: " << stateStream.str() << "\n";
  }
  */

  return {logLike, bgLogLike};
}

double AlignmentModel::logLikelihood(const ReadPair& hit, Transcript& ref) {
  double logLike = salmon::math::LOG_1;
  double bg = salmon::math::LOG_1;
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return logLike;
  }

  if (!hit.isPaired()) {
    if (hit.isLeftOrphan()) {
      auto alnLogProb = logLikelihood(hit.read1, ref, transitionProbsLeft_);
      return alnLogProb.fg - alnLogProb.bg;
    } else {
      auto alnLogProb = logLikelihood(hit.read1, ref, transitionProbsRight_);
      return alnLogProb.fg - alnLogProb.bg;
    }
  }

  bam_seq_t* leftRead =
      (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
  bam_seq_t* rightRead =
      (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read2 : hit.read1;

  size_t leftLen = static_cast<size_t>(bam_seq_len(leftRead));
  size_t rightLen = static_cast<size_t>(bam_seq_len(rightRead));

  if (leftRead) {
    auto alnLogProb = logLikelihood(leftRead, ref, transitionProbsLeft_);
    logLike += alnLogProb.fg;
    bg += alnLogProb.bg;
  }

  if (rightRead) {
    auto alnLogProb = logLikelihood(rightRead, ref, transitionProbsRight_);
    logLike += alnLogProb.fg;
    bg += alnLogProb.bg;
  }

  if (logLike == salmon::math::LOG_0) {
    std::lock_guard<std::mutex> lock(outputMutex_);
    std::cerr << "orphan status: " << hit.orphanStatus << "\n";
    std::cerr << "error likelihood: " << logLike << "\n";
  }

  return logLike - bg;
}

double AlignmentModel::logLikelihood(const UnpairedRead& hit, Transcript& ref) {
  double logLike = salmon::math::LOG_1;
  double bg = salmon::math::LOG_1;
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return logLike;
  }

  bam_seq_t* read = hit.read;
  size_t readLen = static_cast<size_t>(bam_seq_len(read));
  // NOTE: Raise a warning in this case?
  /*
  if (BOOST_UNLIKELY(readLen > maxExpectedLen_)) {
      return logLike;
  }
  */
  auto alnLogProb = logLikelihood(read, ref, transitionProbsLeft_);
  logLike += alnLogProb.fg;
  bg += alnLogProb.bg;

  if (logLike == salmon::math::LOG_0) {
    std::lock_guard<std::mutex> lock(outputMutex_);
    std::cerr << "error log likelihood: " << logLike << "\n";
  }

  return logLike - bg;
}

void AlignmentModel::update(const UnpairedRead& hit, Transcript& ref, double p,
                            double mass) {
  if (mass == salmon::math::LOG_0) {
    return;
  }
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return;
  }
  bam_seq_t* leftRead = hit.read;
  update(leftRead, ref, p, mass, transitionProbsLeft_);
}

void AlignmentModel::update(
    bam_seq_t* read, Transcript& ref, double p, double mass,
    std::vector<AtomicMatrix<double>>& transitionProbs) {
  using namespace salmon::stringtools;
  bool useQual{false};
  int32_t readIdx{0};
  auto transcriptIdx = bam_pos(read);
  size_t transcriptLen = ref.RefLength;

  // if the read starts before the beginning of the transcript,
  // only consider the part overlapping the transcript
  if (transcriptIdx < 0) {
    readIdx = -transcriptIdx;
    transcriptIdx = 0;
  }
  // unsigned version of transcriptIdx
  size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

  uint32_t* cigar = bam_cigar(read);
  uint32_t cigarLen = bam_cigar_len(read);
  uint8_t* qseq = reinterpret_cast<uint8_t*>(bam_seq(read));
  uint8_t* qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));
  int32_t readLen = bam_seq_len(read);

  if (cigarLen > 0 and cigar) {

    salmon::stringtools::strand readStrand =
        salmon::stringtools::strand::forward;

    bool advanceInRead{false};
    bool advanceInReference{false};
    uint32_t readPosBin{0};
    uint32_t cigarIdx{0};
    uint32_t prevStateIdx{startStateIdx};
    uint32_t curStateIdx{0};
    double invLen = static_cast<double>(readBins_) / readLen;

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
              logger_->warn("(in update()) CIGAR string for read [{}] "
                            "seems inconsistent. It refers to non-existant "
                            "positions in the read!",
                            bam_name(read));
              std::stringstream cigarStream;
              for (size_t j = 0; j < cigarLen; ++j) {
                uint32_t opLen = cigar[j] >> BAM_CIGAR_SHIFT;
                enum cigar_op op =
                    static_cast<enum cigar_op>(cigar[j] & BAM_CIGAR_MASK);
                cigarStream << opLen << opToChr(op);
              }
              logger_->warn("CIGAR = {}", cigarStream.str());
            }
            return;
          }

          curReadBase = samToTwoBit[bam_seqi(qseq, readIdx)];
          readPosBin = static_cast<uint32_t>((readIdx * invLen));
          advanceInRead = false;
        }
        if (advanceInReference) {
          // Shouldn't happen!
          if (uTranscriptIdx >= transcriptLen) {
            if (logger_) {
              logger_->warn(
                  "(in update()) CIGAR string for read [{}] "
                  "seems inconsistent. It refers to non-existant "
                  "positions in the reference! Transcript name "
                  "is {}, length is {}, id is {}. Read things refid is {}",
                  bam_name(read), ref.RefName, transcriptLen, ref.id,
                  bam_ref(read));
            }
            return;
          }

          curRefBase = samToTwoBit[ref.baseAt(uTranscriptIdx, readStrand)];
          advanceInReference = false;
        }

        setBasesFromCIGAROp_(
            op, curRefBase,
            curReadBase); //, readStream, matchStream, refStream);
        curStateIdx = curRefBase * numStates + curReadBase;

        // update the state in the actual model
        transitionProbs[readPosBin].increment(prevStateIdx, curStateIdx,
                                              mass + p);

        prevStateIdx = curStateIdx;
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
  } // if we had a cigar string
}

void AlignmentModel::update(const ReadPair& hit, Transcript& ref, double p,
                            double mass) {
  if (mass == salmon::math::LOG_0) {
    return;
  }
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return;
  }

  if (hit.isPaired()) {
    bam_seq_t* leftRead =
        (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read1 : hit.read2;
    bam_seq_t* rightRead =
        (bam_pos(hit.read1) < bam_pos(hit.read2)) ? hit.read2 : hit.read1;
    update(leftRead, ref, p, mass, transitionProbsLeft_);
    update(rightRead, ref, p, mass, transitionProbsRight_);
  } else if (hit.isLeftOrphan()) {
    bam_seq_t* read = hit.read1;
    update(read, ref, p, mass, transitionProbsLeft_);
  } else if (hit.isRightOrphan()) {
    bam_seq_t* read = hit.read1;
    update(read, ref, p, mass, transitionProbsRight_);
  }
}

// CIGAR string with printing
/*
    std::stringstream readStr;
    std::stringstream matchStr;
    std::stringstream refStr;
    std::stringstream stateStr;

switch (op) {
                    case BAM_UNKNOWN:
                        std::cerr << "ENCOUNTERED UNKNOWN SYMBOL IN CIGAR
STRING!\n"; break; case BAM_CMATCH: readStr << twoBitToChar[curReadBase];
                        matchStr << ((curReadBase == curRefBase) ? ' ' : 'X');
                        refStr << twoBitToChar[curRefBase];
                        // do nothing
                        break;
                    case BAM_CBASE_MATCH:
                        readStr << twoBitToChar[curReadBase];
                        matchStr << ' ';
                        refStr << twoBitToChar[curRefBase];
                        // do nothing
                        break;
                    case BAM_CBASE_MISMATCH:
                        readStr << twoBitToChar[curReadBase];
                        matchStr << 'X';
                        refStr << twoBitToChar[curRefBase];
                        // do nothing
                        break;
                    case BAM_CINS:
                        readStr << twoBitToChar[curReadBase];
                        matchStr << ' ';
                        refStr << '-';
                        curRefBase = ALN_DASH;
                        break;
                    case BAM_CDEL:
                        readStr << '-';
                        matchStr << ' ';
                        refStr << twoBitToChar[curRefBase];
                        curReadBase = ALN_DASH;
                        break;
                    case BAM_CREF_SKIP:
                        readStr << 'N';
                        matchStr << ' ';
                        refStr << twoBitToChar[curRefBase];
                        curReadBase = ALN_REF_SKIP;
                        break;
                    case BAM_CSOFT_CLIP:
                        readStr << twoBitToChar[curReadBase];
                        matchStr << ' ';
                        refStr << 'S';
                        curRefBase = ALN_SOFT_CLIP;
                        break;
                    case BAM_CHARD_CLIP:
                        readStr << 'H';
                        matchStr << ' ';
                        refStr << 'H';
                        curRefBase = ALN_HARD_CLIP;
                        curReadBase = ALN_HARD_CLIP;
                        break;
                    case BAM_CPAD:
                        readStr << 'P';
                        matchStr << ' ';
                        refStr << 'P';
                        curRefBase = ALN_PAD;
                        curReadBase = ALN_PAD;
                        break;
                }

                curStateIdx = curRefBase * numStates + curReadBase;


               stateStr << curStateIdx << ", ";
        {
            std::lock_guard<std::mutex> l(outputMutex_);
            std::cerr << "\n\nread:   " << readStr.str() << "\n";
            std::cerr << "        " << matchStr.str() << "\n";
            std::cerr << "ref:    " << refStr.str() << "\n";
            std::cerr << "states: " << stateStr.str() << "\n";
        }
*/
