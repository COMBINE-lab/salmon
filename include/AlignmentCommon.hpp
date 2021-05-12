#ifndef __ALIGNMENTCOMMON_H__
#define __ALIGNMENTCOMMON_H__

#include <atomic>
#include <memory>
#include <mutex>


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

  template<typename T>
  static int32_t alnLen(const T& aln, const T& primary) {
    const auto l = aln.readLen();
    return l != 0 ? l : primary.readLen();
  }

  struct ErrorCount {
  protected:
    int32_t insertions_, deletions_, mismatches_, matches_;
    int32_t sclips_, hclips_; // soft and hard clips
    int32_t fclips_, bclips_; // clips at front and at back
    friend AlignmentCommon;

  public:
    inline int32_t clips() const { return sclips_ + hclips_; }
    // Indels + mismatches
    inline int32_t ims() const { return insertions_ + deletions_ + mismatches_; }
    // Should be equal to the length of the query sequence
    inline int32_t length() const { return insertions_ + mismatches_ + matches_ + sclips_; }
    void clear() {
      insertions_ = deletions_ = mismatches_ = matches_ = sclips_ = hclips_ = fclips_ = bclips_ = 0;
    }
    inline int32_t insertions() const { return insertions_; }
    inline int32_t deletions() const { return deletions_; }
    inline int32_t mismatches() const { return mismatches_; }
    inline int32_t matches() const { return matches_; }
    inline int32_t sclips() const { return sclips_; }
    inline int32_t hclips() const { return hclips_; }
    inline int32_t fclips() const { return fclips_; }
    inline int32_t bclips() const { return bclips_; }

    bool computeErrorCount(bam_seq_t* read, bam_seq_t* primary, Transcript& ref,
                           const char* src);
  };
  bool computeErrorCount(bam_seq_t* read, bam_seq_t* primary, Transcript& ref,
                         ErrorCount& counts, const char* src);

  std::shared_ptr<spdlog::logger> logger_;
  std::atomic<bool> burnedIn_;

  std::mutex throwMutex_;
};

#endif /* __ALIGNMENTCOMMON_H__ */
