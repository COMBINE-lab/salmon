#ifndef ALIGNMENT_MODEL
#define ALIGNMENT_MODEL

#include <mutex>

#include "AlignmentCommon.hpp"
#include "AtomicMatrix.hpp"
#include "oneapi/tbb/concurrent_vector.h"

class AlignmentModel
  : public AlignmentCommon
{
public:
  AlignmentModel(double alpha, uint32_t readBins = 4);

  /**
   *  For unpaired reads, update the error model to account for errors
   *  we've observed in this read pair. primaryAln contains the first
   *  alignment in the alignment group.
   */
  void update(const UnpairedRead& aln, const UnpairedRead& primaryAln,
              Transcript& ref, double p, double mass);

  /**
   * Compute the log-likelihood of the observed unpaired alignment
   * given the current error model. primaryAln contains the first
   * alignment in the alignment group.
   */
  double logLikelihood(const UnpairedRead&, const UnpairedRead& primaryAln, Transcript& ref);

  // The primaryAln is ignored by the ReadPair overlaods (at least for now)

  /**
   *  For paired-end reads, update the error model to account for
   *  errors we've observed in this read pair.
   */
  void update(const ReadPair& aln, const ReadPair& primaryAln,
              Transcript& ref, double p, double mass);

  /**
   * Compute the log-likelihood of the observed paire-end alignment given the
   * current error model.
   */
  double logLikelihood(const ReadPair& aln, const ReadPair& primaryAln, Transcript& ref);

  void normalize();

private:

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
  void update(bam_seq_t* read, bam_seq_t* primary, Transcript& ref, double p, double mass,
              std::vector<AtomicMatrix<double>>& mismatchProfile);
  AlnModelProb logLikelihood(bam_seq_t* read, bam_seq_t* primary, Transcript& ref,
                       std::vector<AtomicMatrix<double>>& mismatchProfile);

  // NOTE: Do these need to be concurrent_vectors as before?
  // Store the mismatch probability tables for the left and right reads
  std::vector<AtomicMatrix<double>> transitionProbsLeft_;
  std::vector<AtomicMatrix<double>> transitionProbsRight_;

  bool isEnabled_;
  // size_t maxLen_;
  size_t readBins_;
  // Maintain a mutex in case the error model wants to talk to the
  // console / log.
  std::mutex outputMutex_;
};

#endif // ERROR_MODEL
