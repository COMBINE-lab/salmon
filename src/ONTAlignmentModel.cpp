#include <cstdio>
#include <iostream>
#include <map>
#include <tuple>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY

#include "ONTAlignmentModel.hpp"
#include "SalmonMath.hpp"
#include "SalmonStringUtils.hpp"
#include "Transcript.hpp"
#include "UnpairedRead.hpp"

ONTAlignmentModel::ONTAlignmentModel(double alpha, uint32_t readBins)
    : isEnabled_(true)
    , readBins_(readBins)
    , errorModel_(maxReadLen / binLen + 1)
{
  for(auto& avg: errorModel_) {
    avg.number = 0;
    avg.sum    = 0.0;
  }
}

// ONTAlignmentModel::AlnModelProb ONTAlignmentModel::logLikelihood(bam_seq_t* read, Transcript& ref,
//                                                                  std::vector<AtomicMatrix<double>>& transitionProbs)
// {
//   throw std::runtime_error("ONTAlignmentModel not yet implemented");
//   return {0.0, 0.0};
// }

double ONTAlignmentModel::logLikelihood(const UnpairedRead& hit, Transcript& ref) {
  std::cerr << "ONTAlignmentModel log likelihood not yet implemented" << std::endl;
  throw std::runtime_error("ONTAlignmentModel log likelihood not yet implemented");
  double logLike = salmon::math::LOG_1;
  double bg = salmon::math::LOG_1;
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return logLike;
  }
  return logLike - bg;
}

// Update the probability model. The reads are binned based on their
// length (the length after soft clipping). For each bin, we assume a
// Binomial distribution B(p,n) (where n is the length of the read).
//
// The probability p in each bin is estimated as the empirical mean of
// the number of errors in the reads (all reads counted the same at
// this point: indels and mutations).
void ONTAlignmentModel::update(const UnpairedRead& hit, Transcript& ref, double p,
                               double mass) {
  if (mass == salmon::math::LOG_0) {
    return;
  }
  if (BOOST_UNLIKELY(!isEnabled_)) {
    return;
  }

  using namespace salmon::stringtools;
  bool useQual{false};
  int32_t readIdx{0};
  auto         transcriptIdx = bam_pos(hit.read);
  const size_t transcriptLen = ref.RefLength;

  // if the read starts before the beginning of the transcript,
  // only consider the part overlapping the transcript
  if (transcriptIdx < 0) {
    readIdx = -transcriptIdx;
    transcriptIdx = 0;
  }
  // unsigned version of transcriptIdx
  size_t uTranscriptIdx = static_cast<size_t>(transcriptIdx);

  uint32_t*      cigar       = bam_cigar(hit.read);
  const uint32_t cigarLen    = bam_cigar_len(hit.read);
  uint8_t*       qseq        = reinterpret_cast<uint8_t*>(bam_seq(hit.read));
  //  uint8_t*       qualStr = reinterpret_cast<uint8_t*>(bam_qual(read));
  const int32_t  readLen     = bam_seq_len(hit.read);

  if (cigarLen == 0 || !cigar)
    return;

  const salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
  bool advanceInRead{false};
  bool advanceInReference{false};

  uint32_t errors = 0;
  uint32_t clips = 0;

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
                          bam_name(hit.read));
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
                          bam_name(hit.read), ref.RefName, transcriptLen, ref.id,
                          bam_ref(hit.read));
          }
          return;
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

  int32_t alignLen = readLen - clips;
  double errorRate = (double)errors / alignLen;
  if(errorRate > 1.0) { // Should not happen
    if (logger_) {
      logger_->warn("(in update()) CIGAR string for read [{}] "
                    "seems inconsistent. It implied an error rate "
                    "greater than 1",
                    bam_name(hit.read));
    }
    return;
  }
  int32_t bin = std::min(alignLen / binLen, (uint32_t)errorModel_.size() - 1);
  ++errorModel_[bin].number;
  salmon::utils::incLoop(errorModel_[bin].sum, errorRate);
}

// void ONTAlignmentModel::update(
//     bam_seq_t* read, Transcript& ref, double p, double mass,
//     std::vector<AtomicMatrix<double>>& transitionProbs) {
//   throw std::runtime_error("ONTAlignmentModel not yet implemented");
// }

void ONTAlignmentModel::print_model(std::ostream& os) {
  for(size_t i = 0; i < errorModel_.size(); ++i) {
    os << (i * binLen) << " - " << ((i+1) * binLen) << ' ' << (errorModel_[i].sum / errorModel_[i].number) << '\n';
  }
}
