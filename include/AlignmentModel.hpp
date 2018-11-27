#ifndef ALIGNMENT_MODEL
#define ALIGNMENT_MODEL

#include <atomic>
#include <memory>
#include <mutex>

// logger includes
#include "spdlog/spdlog.h"

#include "AtomicMatrix.hpp"
#include "tbb/concurrent_vector.h"

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

struct UnpairedRead;
struct ReadPair;
class Transcript;

class AlignmentModel {
public:
  AlignmentModel(double alpha, uint32_t readBins = 4);
  bool burnedIn();
  void burnedIn(bool burnedIn);

  /**
   *  For unpaired reads, update the error model to account
   *  for errors we've observed in this read pair.
   */
  void update(const UnpairedRead&, Transcript& ref, double p, double mass);

  /**
   * Compute the log-likelihood of the observed unpaired alignment given the
   * current error model.
   */
  double logLikelihood(const UnpairedRead&, Transcript& ref);

  /**
   *  For paired-end reads, update the error model to account
   *  for errors we've observed in this read pair.
   */
  void update(const ReadPair&, Transcript& ref, double p, double mass);

  /**
   * Compute the log-likelihood of the observed paire-end alignment given the
   * current error model.
   */
  double logLikelihood(const ReadPair&, Transcript& ref);

  bool hasIndel(UnpairedRead& r);
  bool hasIndel(ReadPair& r);

  void setLogger(std::shared_ptr<spdlog::logger> logger);
  bool hasLogger();

  void normalize();

private:
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

  void setBasesFromCIGAROp_(enum cigar_op op, size_t& curRefBase,
                            size_t& curReadBase);
  // std::stringstream& readStr, std::stringstream& matchStr,
  // std::stringstream& refstr);

  // 11 states --- the 9 listed in "AlignmentModelChar" plus START (9) i
  // and END (10)
  constexpr static uint32_t numAlignmentStates() { return 82; }
  constexpr static uint32_t numStates = 9;
  constexpr static uint32_t startStateIdx = 81;

  struct AlnModelProb {
    double fg{salmon::math::LOG_1};
    double bg{salmon::math::LOG_1};
  };

  /**
   * These functions, which work directly on bam_seq_t* types, drive the
   * update() and logLikelihood() methods above.
   */
  void update(bam_seq_t* read, Transcript& ref, double p, double mass,
              std::vector<AtomicMatrix<double>>& mismatchProfile);
  AlnModelProb logLikelihood(bam_seq_t* read, Transcript& ref,
                       std::vector<AtomicMatrix<double>>& mismatchProfile);
  bool hasIndel(bam_seq_t* r);

  // NOTE: Do these need to be concurrent_vectors as before?
  // Store the mismatch probability tables for the left and right reads
  std::vector<AtomicMatrix<double>> transitionProbsLeft_;
  std::vector<AtomicMatrix<double>> transitionProbsRight_;

  bool isEnabled_;
  // size_t maxLen_;
  size_t readBins_;
  std::atomic<bool> burnedIn_;
  std::shared_ptr<spdlog::logger> logger_;
  // Maintain a mutex in case the error model wants to talk to the
  // console / log.
  std::mutex outputMutex_;
};

#endif // ERROR_MODEL
