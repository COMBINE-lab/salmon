#ifndef ONT_ALIGNMENT_MODEL
#define ONT_ALIGNMENT_MODEL

#include <atomic>
#include <memory>
#include <mutex>

// logger includes
#include "spdlog/spdlog.h"

#include "AlignmentCommon.hpp"

#include "AtomicMatrix.hpp"
// #include "tbb/concurrent_vector.h"


class ONTAlignmentModel
  : public AlignmentCommon
{
public:
  static const uint32_t maxReadLen = 50000; // XXX: That should be a paramater. Read longer than that are binned together
  static const uint32_t binLen = 100; // XXX: That should be a parameter
  
   /**
    * kmerLength specifies the length of the kmers that each state represents in the Markov Model
    * stepSize specifies how far we move the kmer window when looking at state transition 
    * (ie stepSize = kmerLength means no overlap of kmers,stepSize = kmerLength-1 would mean kmers that overlap by 1 etc. ) 
    * We must have 0 < stepSize <= kmerLength.
    */
  static const uint32_t kmerLength = 50;
  static const uint32_t stepSize = 50;
 
  ONTAlignmentModel(double alpha, uint32_t readBins = 4);
  ~ONTAlignmentModel() {   }

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
  //added for new ONT model
   struct AlnModelProb {
    double fg{salmon::math::LOG_1};
    double bg{salmon::math::LOG_1};
  };
  AlnModelProb logLikelihood(bam_seq_t* read, bam_seq_t* primary, Transcript& ref,
                       std::vector<AtomicMatrix<double>>& mismatchProfile);
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
  
  void update(bam_seq_t* read, bam_seq_t* primary, Transcript& ref, double p, double mass,
            std::vector<AtomicMatrix<double>>& mismatchProfile);

  //Transtion probabilities for the mismatch/indel error Markov Model
  std::vector<AtomicMatrix<double>> transitionProbs_;
};

#endif // ERROR_MODEL
