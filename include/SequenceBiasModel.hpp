#ifndef SEQUENCE_BIAS_MODEL
#define SEQUENCE_BIAS_MODEL

#include <atomic>
#include <memory>
#include <mutex>

// logger includes
#include "spdlog/spdlog.h"

#include "AtomicMatrix.hpp"
#include "tbb/concurrent_vector.h"

class Transcript;
class LibraryFormat;

class SequenceBiasModel {
public:
  SequenceBiasModel(double alpha,
                    uint32_t windowSize = 15 // must be odd
  );

  bool burnedIn();
  void burnedIn(bool burnedIn);

  /**
   *  For unpaired reads, update the bias model to account
   *  for this read
   */
  bool update(Transcript& ref, int32_t pos, LibraryFormat libFormat,
              double mass, double prob);

  /**
   *  For unpaired reads, update the bias model to account
   *  for this read
   */
  bool update(Transcript& ref, int32_t pos, int32_t pos2,
              LibraryFormat libFormat, double mass, double prob);

  /**
   *
   *
   */
  double biasFactor(Transcript& ref, int32_t pos, LibraryFormat libFmt);

  /**
   *
   *
   */
  double biasFactor(Transcript& ref, int32_t pos1, int32_t pos2,
                    LibraryFormat libFormat);

  void setLogger(std::shared_ptr<spdlog::logger> logger);
  bool hasLogger();
  std::string toString();
  // void normalize();

private:
  double seqProb(Transcript& ref, int32_t pos, bool isFwd,
                 AtomicMatrix<double>& profile);
  /**
   * These functions, which work directly on bam_seq_t* types, drive the
   * update() and logLikelihood() methods above.
   */
  bool update(Transcript& ref, int32_t pos, bool fwd, double mass, double prob,
              AtomicMatrix<double>& sequenceModel);

  double biasFactor(Transcript& ref, int32_t pos, LibraryFormat libFormat,
                    AtomicMatrix<double>& foregroundModel,
                    AtomicMatrix<double>& backgroundModel);

  // A, C, G, T
  constexpr static uint32_t numBases() { return 4; }

  // NOTE: Do these need to be concurrent_vectors as before?
  // Store the mismatch probability tables for the left and right reads
  AtomicMatrix<double> biasLeftForeground_;
  AtomicMatrix<double> biasLeftBackground_;
  AtomicMatrix<double> biasRightForeground_;
  AtomicMatrix<double> biasRightBackground_;

  uint32_t windowSize_;
  //bool isEnabled_;
  std::atomic<bool> burnedIn_;
  std::shared_ptr<spdlog::logger> logger_;
  // Maintain a mutex in case we want to talk to the
  // console / log.
  std::mutex outputMutex_;
};

#endif // SEQUENCE_BIAS_MODEL
