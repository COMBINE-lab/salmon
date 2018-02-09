#include <cstdio>
#include <iostream>

#include <boost/config.hpp> // for BOOST_LIKELY/BOOST_UNLIKELY

#include "LibraryFormat.hpp"
#include "SequenceBiasModel.hpp"
#include "Transcript.hpp"
#include "spdlog/fmt/fmt.h"
#include "spdlog/fmt/ostr.h"

SequenceBiasModel::SequenceBiasModel(double alpha, uint32_t windowSize)
    : // Let's try a 0th order model first
      biasLeftForeground_(AtomicMatrix<double>(windowSize, numBases(), alpha)),
      biasLeftBackground_(AtomicMatrix<double>(windowSize, numBases(), alpha)),
      biasRightForeground_(AtomicMatrix<double>(windowSize, numBases(), alpha)),
      biasRightBackground_(AtomicMatrix<double>(windowSize, numBases(), alpha)),
      //isEnabled_(true),
      windowSize_(windowSize), burnedIn_(false) {}

bool SequenceBiasModel::burnedIn() { return burnedIn_; }
void SequenceBiasModel::burnedIn(bool burnedIn) { burnedIn_ = burnedIn; }

void SequenceBiasModel::setLogger(std::shared_ptr<spdlog::logger> logger) {
  logger_ = logger;
}

bool SequenceBiasModel::hasLogger() { return (logger_) ? true : false; }

bool SequenceBiasModel::update(Transcript& ref, int32_t pos, bool fwd,
                               double mass, double prob,
                               AtomicMatrix<double>& sequenceProfile) {

  int32_t halfWindow = windowSize_ >> 1;
  int32_t txpLen = ref.RefLength;

  // We must be at least halfWindow bases into the transcript
  // and no fewer than halfWindow bases from the end
  if (pos < halfWindow) {
    return false;
  }
  if (pos + halfWindow >= txpLen) {
    return false;
  }

  int32_t start, stop, inc;
  salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
  if (fwd) {
    start = pos - halfWindow;
    stop = pos + halfWindow;
    for (int p = start; p <= stop; ++p) {
      size_t currRefBase =
          salmon::stringtools::samToTwoBit[ref.baseAt(p, readStrand)];
      sequenceProfile.increment(p - start, currRefBase, mass + prob);
    }
  } else {
    start = pos + halfWindow;
    stop = pos - halfWindow;
    for (int p = start; p >= stop; --p) {
      size_t currRefBase =
          salmon::stringtools::samToTwoBit[ref.baseAt(p, readStrand)];
      sequenceProfile.increment(start - p, currRefBase, mass + prob);
    }
  }
  return true;
}

bool SequenceBiasModel::update(Transcript& ref, int32_t pos,
                               LibraryFormat libFmt, double mass, double prob) {

  uint32_t offset{100};
  int32_t bgPos = pos - offset;
  bool fwd = true;
  if (libFmt.strandedness == ReadStrandedness::S) {
    bgPos = pos - offset;
    fwd = true;
  } else {
    bgPos = pos + offset;
    fwd = false;
  }
  int32_t halfWindow = windowSize_ >> 1;
  int32_t txpLen = ref.RefLength;

  // We must be at least halfWindow bases into the transcript
  // and no fewer than halfWindow bases from the end
  if (pos - halfWindow < 100) {
    return false;
  }
  if (pos + halfWindow >= txpLen - 100) {
    return false;
  }
  if (bgPos - halfWindow < 100) {
    return false;
  }
  if (bgPos + halfWindow >= txpLen - 100) {
    return false;
  }

  update(ref, pos, fwd, mass, prob, biasLeftForeground_);
  update(ref, bgPos, fwd, mass, prob, biasLeftBackground_);
  return true;
}

bool SequenceBiasModel::update(Transcript& ref, int32_t pos1, int32_t pos2,
                               LibraryFormat libFmt, double mass, double prob) {

  int32_t leftPos, rightPos;
  if (pos1 < pos2) {
    leftPos = pos1;
    rightPos = pos2;
  } else {
    rightPos = pos1;
    leftPos = pos2;
  }

  uint32_t offset{100};
  int32_t bgPosLeft = leftPos - offset;
  int32_t bgPosRight = rightPos + offset;
  int32_t halfWindow = windowSize_ >> 1;
  int32_t txpLen = ref.RefLength;

  // We must be at least halfWindow bases into the transcript
  // and no fewer than halfWindow bases from the end
  if (leftPos < halfWindow) {
    return false;
  }
  if (leftPos + halfWindow >= txpLen) {
    return false;
  }
  if (bgPosLeft < halfWindow) {
    return false;
  }
  if (bgPosLeft + halfWindow >= txpLen) {
    return false;
  }
  if (rightPos < halfWindow) {
    return false;
  }
  if (rightPos + halfWindow >= txpLen) {
    return false;
  }
  if (bgPosRight < halfWindow) {
    return false;
  }
  if (bgPosRight + halfWindow >= txpLen) {
    return false;
  }

  update(ref, leftPos, true, mass, prob, biasLeftForeground_);
  update(ref, bgPosLeft, true, mass, prob, biasLeftBackground_);
  update(ref, rightPos, true, mass, prob, biasRightForeground_);
  update(ref, bgPosRight, true, mass, prob, biasRightBackground_);
  return true;
}

double SequenceBiasModel::seqProb(Transcript& ref, int32_t pos, bool isFwd,
                                  AtomicMatrix<double>& profile) {
  int32_t halfWindow = windowSize_ >> 1;
  salmon::stringtools::strand readStrand = salmon::stringtools::strand::forward;
  int32_t start, stop;
  double prob = salmon::math::LOG_1;

  if (isFwd) {
    start = pos - halfWindow;
    stop = pos + halfWindow;
    for (int p = start; p <= stop; ++p) {
      size_t currRefBase =
          salmon::stringtools::samToTwoBit[ref.baseAt(p, readStrand)];
      prob += profile(p - start, currRefBase);
    }
  } else {
    start = pos + halfWindow;
    stop = pos - halfWindow;
    for (int p = start; p >= stop; --p) {
      size_t currRefBase =
          salmon::stringtools::samToTwoBit[ref.baseAt(p, readStrand)];
      prob += profile(start - p, currRefBase);
    }
  }

  return prob;
}

double SequenceBiasModel::biasFactor(Transcript& ref, int32_t pos,
                                     LibraryFormat libFormat,
                                     AtomicMatrix<double>& foregroundProfile,
                                     AtomicMatrix<double>& backgroundProfile) {

  int32_t halfWindow = windowSize_ >> 1;
  int32_t txpLen = ref.RefLength;

  bool isFwd = true;
  if (libFormat.strandedness == ReadStrandedness::S) {
    isFwd = true;
  } else {
    isFwd = false;
  }

  // We must be at least halfWindow bases into the transcript
  // and no fewer than halfWindow bases from the end
  if (pos - halfWindow > 100) {
    return salmon::math::LOG_1;
  }
  if (pos + halfWindow >= txpLen - 100) {
    return salmon::math::LOG_1;
  }

  double condProbBack = seqProb(ref, pos, isFwd, backgroundProfile);
  double condProbFore = seqProb(ref, pos, isFwd, foregroundProfile);
  return (condProbBack - condProbFore);
}

double SequenceBiasModel::biasFactor(Transcript& ref, int32_t pos1,
                                     int32_t pos2, LibraryFormat libFormat) {

  int32_t leftPos, rightPos;
  if (pos1 < pos2) {
    leftPos = pos1;
    rightPos = pos2;
  } else {
    rightPos = pos1;
    leftPos = pos2;
  }

  int32_t halfWindow = windowSize_ >> 1;
  int32_t txpLen = ref.RefLength;

  // We must be at least halfWindow bases into the transcript
  // and no fewer than halfWindow bases from the end
  if (leftPos < halfWindow) {
    return salmon::math::LOG_1;
  }
  if (leftPos + halfWindow >= txpLen) {
    return salmon::math::LOG_1;
  }
  if (rightPos < halfWindow) {
    return salmon::math::LOG_1;
  }
  if (rightPos + halfWindow >= txpLen) {
    return salmon::math::LOG_1;
  }

  double biasLeft = biasFactor(ref, leftPos, libFormat, biasLeftForeground_,
                               biasLeftBackground_);
  double biasRight = biasFactor(ref, rightPos, libFormat, biasRightForeground_,
                                biasRightBackground_);

  // This is in log space, so it's actually log(leftBias * rightBias).  Is
  // multiplying the biases the right thing to do here, or should we average
  // them (i.e. log([leftBias + rightBias]/2))?
  return biasLeft + biasRight;
}

double SequenceBiasModel::biasFactor(Transcript& ref, int32_t pos,
                                     LibraryFormat libFormat) {

  int32_t halfWindow = windowSize_ >> 1;
  int32_t txpLen = ref.RefLength;

  // We must be at least halfWindow bases into the transcript
  // and no fewer than halfWindow bases from the end
  if (pos < halfWindow) {
    return salmon::math::LOG_1;
  }
  if (pos + halfWindow >= txpLen) {
    return salmon::math::LOG_1;
  }

  double bias =
      biasFactor(ref, pos, libFormat, biasLeftForeground_, biasLeftBackground_);
  return bias;
}

std::string SequenceBiasModel::toString() {
  fmt::MemoryWriter srep;
  srep << "{";
  srep << "\"fg\": {";
  for (size_t i = 0; i < windowSize_; ++i) {
    srep << "\"" << i << "\":{";
    for (size_t r = 0; r < numBases(); ++r) {
      switch (r) {
      case 0:
        srep << "\"A\":";
        break;
      case 1:
        srep << "\"C\":";
        break;
      case 2:
        srep << "\"G\":";
        break;
      case 3:
        srep << "\"T\":";
        break;
      }
      srep << std::exp(biasLeftForeground_(i, r));
      if (r < numBases() - 1) {
        srep << ", ";
      } else {
        srep << "} ";
      }
    }
    if (i < windowSize_ - 1) {
      srep << ",";
    }
  }

  srep << "}, ";
  srep << "\"bg\": {";
  for (size_t i = 0; i < windowSize_; ++i) {
    srep << "\"" << i << "\":{";
    for (size_t r = 0; r < numBases(); ++r) {
      switch (r) {
      case 0:
        srep << "\"A\":";
        break;
      case 1:
        srep << "\"C\":";
        break;
      case 2:
        srep << "\"G\":";
        break;
      case 3:
        srep << "\"T\":";
        break;
      }
      srep << std::exp(biasLeftBackground_(i, r));
      if (r < numBases() - 1) {
        srep << ", ";
      } else {
        srep << "} ";
      }
    }
    if (i < windowSize_ - 1) {
      srep << ",";
    }
  }
  srep << "}}";

  return srep.str();
}
