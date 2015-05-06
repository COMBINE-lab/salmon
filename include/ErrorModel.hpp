#ifndef ERROR_MODEL
#define ERROR_MODEL

#include <mutex>
#include <atomic>
#include "tbb/concurrent_vector.h"
#include "AtomicMatrix.hpp"

extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
#undef max
#undef min
}

struct UnpairedRead;
struct ReadPair;
class Transcript;

class ErrorModel {
public:
    ErrorModel(double alpha, uint32_t maxExpectedLen);
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

private:
    /**
     * These functions, which work directly on bam_seq_t* types, drive the
     * update() and logLikelihood() methods above.
     */
    void update(bam_seq_t* read, Transcript& ref, double p, double mass, std::vector<AtomicMatrix<double>>& mismatchProfile);
    double logLikelihood(bam_seq_t* read, Transcript& ref, std::vector<AtomicMatrix<double>>& mismatchProfile);

    // NOTE: Do these need to be concurrent_vectors as before?
    // Store the mismatch probability tables for the left and right reads
    std::vector<AtomicMatrix<double>> mismatchLeft_;
    std::vector<AtomicMatrix<double>> mismatchRight_;
    //tbb::concurrent_vector<AtomicMatrix<double>> mismatchLeft_;
    //tbb::concurrent_vector<AtomicMatrix<double>> mismatchRight_;

    bool isEnabled_;
    size_t maxExpectedLen_;
    size_t maxLen_;
    std::atomic<bool> burnedIn_;
    // Maintain a mutex in case the error model wants to talk to the
    // console / log.
    std::mutex outputMutex_;
};

#endif // ERROR_MODEL
