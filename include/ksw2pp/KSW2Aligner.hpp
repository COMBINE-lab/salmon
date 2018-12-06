#ifndef __KSW2_ALIGNER_HPP__
#define __KSW2_ALIGNER_HPP__

#include <memory>
#include <vector>
#include <cstdlib>

extern "C" {
#include "ksw2pp/kalloc.h"
#include "ksw2pp/ksw2.h"
}

namespace ksw2pp {

/**
 * When we use a unique_ptr to hold a kalloc allocator, this is the
 * deleter we use to call the appropriate function to free / destroy the
 *allocator.
 **/
class KallocDeleter {
public:
  void operator()(void* p) { km_destroy(p); }
};

enum class KSW2AlignmentType : uint8_t { GLOBAL = 1, EXTENSION = 2 };

// Just like Int2Type from
// https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Int-To-Type
template <KSW2AlignmentType I> struct EnumToType {
  enum { value = static_cast<uint8_t>(I) };
};

// A structure to hold the relvant parameters for the aligner
struct KSW2Config {
  int8_t gapo = -1;
  int8_t gape = -1;
  int bandwidth = -1;
  int dropoff = -1;
  int flag = 0;
  int alphabetSize = 5;
  int end_bonus = 10;
  KSW2AlignmentType atype = KSW2AlignmentType::EXTENSION;
};

class KSW2Aligner {

public:
  KSW2Aligner(int8_t match = 2, int8_t mismatch = -4);
  KSW2Aligner(std::vector<int8_t> mat);

  int transformSequenceKSW2(const char* const queryOriginal, const int queryLength,
                            std::vector<unsigned char>& queryTransformed);

  int transformSequencesKSW2(const char* const queryOriginal, const int queryLength,
                             const char* const targetOriginal, const int targetLength,
                             std::vector<unsigned char>& queryTransformed,
                             std::vector<unsigned char>& targetTransformed);
  /**
   * Variants of the operator that require both an explicit type tag to
   * determine the type of alignment to perform, as well as a pointer to
   * an `ksw_extz_t` structure that will be used to hold the output.
   */
  int operator()(const char* const queryOriginal, const int queryLength,
                 const char* const targetOriginal, const int targetLength,
                 ksw_extz_t* ez, EnumToType<KSW2AlignmentType::GLOBAL>);

  int operator()(const char* const queryOriginal, const int queryLength,
                 const char* const targetOriginal, const int targetLength,
                 ksw_extz_t* ez, EnumToType<KSW2AlignmentType::EXTENSION>);

  int operator()(const uint8_t* const queryOriginal, const int queryLength,
                 const uint8_t* const targetOriginal, const int targetLength,
                 ksw_extz_t* ez, EnumToType<KSW2AlignmentType::GLOBAL>);

  int operator()(const uint8_t* const queryOriginal, const int queryLength,
                 const uint8_t* const targetOriginal, const int targetLength,
                 ksw_extz_t* ez, EnumToType<KSW2AlignmentType::EXTENSION>);

  /**
   * Variants of the operator that do not require an output
   * `ksw_extz_t*` variable.  They will store the result in this object's
   * `result_` variable, which can then be queried with the `result()` method.
   */
  int operator()(const char* const queryOriginal, const int queryLength,
                 const char* const targetOriginal, const int targetLength,
                 EnumToType<KSW2AlignmentType::GLOBAL>);

  int operator()(const char* const queryOriginal, const int queryLength,
                 const char* const targetOriginal, const int targetLength,
                 EnumToType<KSW2AlignmentType::EXTENSION>);

  int operator()(const uint8_t* const queryOriginal, const int queryLength,
                 const uint8_t* const targetOriginal, const int targetLength,
                 EnumToType<KSW2AlignmentType::GLOBAL>);

  int operator()(const uint8_t* const queryOriginal, const int queryLength,
                 const uint8_t* const targetOriginal, const int targetLength,
                 EnumToType<KSW2AlignmentType::EXTENSION>);

  /**
   * Variants of the operator that require neither an output
   * `ksw_extz_t*` variable or an explicit type tag to select the type of
   *alignment
   * to perform.  These variants will perform whichever type of alignment is
   * currently specified in the `atype` member of this class's `config_`
   * object.
   **/
  int operator()(const char* const queryOriginal, const int queryLength,
                 const char* const targetOriginal, const int targetLength);

  int operator()(const uint8_t* const queryOriginal, const int queryLength,
                 const uint8_t* const targetOriginal, const int targetLength);

  KSW2Config& config() { return config_; }
  const ksw_extz_t& result() { return result_; }
  void freeCIGAR(ksw_extz_t* ez) {
    if (ez->cigar and kalloc_allocator_) {
      kfree(kalloc_allocator_.get(), ez->cigar);
    }
  }

private:
  std::vector<uint8_t> query_;
  std::vector<uint8_t> target_;
  ksw_extz_t result_;
  std::unique_ptr<void, KallocDeleter> kalloc_allocator_{nullptr,
                                                         KallocDeleter()};
  std::vector<int8_t> mat_;
  KSW2Config config_;
  bool haveSSE41{false};
  bool haveSSE2{false};
};
} // namespace ksw2pp
#endif //__KSW2_ALIGNER_HPP__
