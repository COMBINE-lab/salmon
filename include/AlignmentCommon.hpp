#ifndef __ALIGNMENTCOMMON_H__
#define __ALIGNMENTCOMMON_H__

#include <atomic>
#include <memory>

// logger includes
#include "spdlog/spdlog.h"

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

struct UnpairedRead;
struct ReadPair;
class Transcript;


// Common functionalities to different alignment models
class AlignmentCommon {
public:
  AlignmentCommon()
    : burnedIn_(false)
  { }

  bool burnedIn() { return burnedIn_; }
  void burnedIn(bool burnedIn) { burnedIn_ = burnedIn; }

  void setLogger(std::shared_ptr<spdlog::logger> logger) { logger_ = logger; }
  bool hasLogger() { return (logger_) ? true : false; }


  static bool hasIndel(ReadPair& hit);
  static bool hasIndel(UnpairedRead& hit);

protected:
  enum AlignmentModelChar {
    ALN_A = 0,
    ALN_C = 1,
    ALN_G = 2,
    ALN_T = 3,
    ALN_DASH = 4,
    ALN_SOFT_CLIP = 5,
    ALN_HARD_CLIP = 6,
    ALN_PAD = 7,
    ALN_REF_SKIP = 8
  };

  static bool hasIndel(bam_seq_t* read);
  static void setBasesFromCIGAROp_(enum cigar_op op, size_t& curRefBase, size_t& curReadBase);
  static char opToChr(enum cigar_op op);

  std::shared_ptr<spdlog::logger> logger_;
  std::atomic<bool> burnedIn_;
};

#endif /* __ALIGNMENTCOMMON_H__ */
