#ifndef ONT_ALIGNMENT_MODEL
#define ONT_ALIGNMENT_MODEL

#include <atomic>
#include <memory>
#include <mutex>

// logger includes
#include "spdlog/spdlog.h"

#include "AlignmentCommon.hpp"

// #include "AtomicMatrix.hpp"
// #include "oneapi/tbb/concurrent_vector.h"


class ONTAlignmentModel
  : public AlignmentCommon
{
public:
  static const uint32_t maxReadLen = 50000; // XXX: That should be a paramater. Read longer than that are binned together
  static const uint32_t binLen = 100; // XXX: That should be a parameter

  ONTAlignmentModel(double alpha, uint32_t readBins = 4);
  ~ONTAlignmentModel() { }

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
  double logLikelihood(const UnpairedRead& aln, const UnpairedRead& primaryAln, Transcript& ref);

  void normalize();

  void printModel(std::ostream&);

private:
  // void ONTAlignmentModel::update(bam_seq_t* read, Transcript& ref, double p, double mass,
  //                                std::vector<AtomicMatrix<double>>& transitionProbs);
  bool isEnabled_;
  // size_t maxLen_;
  size_t readBins_;

  // Maintain a mutex in case the error model wants to talk to the
  // console / log.
  bool printed;
  std::mutex outputMutex_;

  struct average {
    std::atomic<double> mass;
    std::atomic<double> sum;
    average() : mass(0.0), sum(0.0) { }
  };
  // Error model. Probability parameter p of the binomial distribution
  // B(p,n) for each read in a bin (based on length n).
  std::vector<average> errorModel_;

  // Clip length model. Geometric distribution with parameter
  // p. Binned for read size.
  // Separate models are considered for front and back clips
  std::vector<average> frontClipModel_;
  std::vector<average> backClipModel_;
};

#endif // ERROR_MODEL
