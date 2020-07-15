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
  static bool hasIndel(bam_seq_t* read);

  std::shared_ptr<spdlog::logger> logger_;
  std::atomic<bool> burnedIn_;
};

#endif /* __ALIGNMENTCOMMON_H__ */
