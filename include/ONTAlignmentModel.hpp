#ifndef ONT_ALIGNMENT_MODEL
#define ONT_ALIGNMENT_MODEL

#include <atomic>
#include <memory>
#include <mutex>

// logger includes
#include "spdlog/spdlog.h"

#include "AlignmentCommon.hpp"

// #include "AtomicMatrix.hpp"
// #include "tbb/concurrent_vector.h"


class ONTAlignmentModel
  : public AlignmentCommon
{
public:
  static const uint32_t maxReadLen = 50000; // XXX: That should be a paramater. Read longer than that are binned together
  static const uint32_t binLen = 100; // XXX: That should be a parameter

  ONTAlignmentModel(double alpha, uint32_t readBins = 4);
  ~ONTAlignmentModel() { }

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

  void normalize();

  void print_model(std::ostream&);

private:
  // void ONTAlignmentModel::update(bam_seq_t* read, Transcript& ref, double p, double mass,
  //                                std::vector<AtomicMatrix<double>>& transitionProbs);
  bool isEnabled_;
  // size_t maxLen_;
  size_t readBins_;
  std::shared_ptr<spdlog::logger> logger_;
  // Maintain a mutex in case the error model wants to talk to the
  // console / log.
  std::mutex outputMutex_;

  // Error model. Probability parameter p of the binomial distribution
  // B(p,n) for each read in a bin (based on length n).
  struct average {
    std::atomic<uint32_t> number;
    std::atomic<double>   sum;
  };
  std::vector<average> errorModel_;
};

#endif // ERROR_MODEL
