#ifndef SALMON_OPTS_HPP
#define SALMON_OPTS_HPP

#include <boost/filesystem.hpp>

// Logger includes
#include "spdlog/spdlog.h"

#include <memory> // for shared_ptr


/**
  * A structure to hold some common options used
  * by Salmon so that we don't have to pass them
  * all around as separate arguments.
  */
struct SalmonOpts {
    // The options below are adopted from the mem_opt_t structure of BWA
    /*
    int maxOccurences; // maximum number of allowable occurences of (S)MEM
    int minSeedOccurences; // try to split a seed into smaller seeds if it occurs
                           // fewer than this many times.
    int minSeedLen; // A seed must be at least this long.
    float splitFactor; // Split a seed if it's longer than splitFactor * minSeedLen.
    int flag; // Used by bwa
    bool maxMEMIntervals; // If true, don't split (S)MEMs into MEMs
    */

    SalmonOpts() : allowOrphans(false), splitSpanningSeeds(false), noFragLengthDist(false),
                   noEffectiveLengthCorrection(false), useReadCompat(false),
                   maxReadOccs(200), extraSeedPass(false),
                   mappingCacheMemoryLimit(5000000) {}

    bool allowOrphans; // Consider orphaned reads when performing lightweight alignemnt.

    bool splitSpanningSeeds; // Attempt to split seeds that span multiple transcripts.

    bool noFragLengthDist; // Do not give a fragment assignment a likelihood based on an emperically
                           // observed fragment length distribution.

    bool noEffectiveLengthCorrection; // Don't take the fragment length distribution into
                                      // account when computing the probability that a
                                     // fragment was generated from a transcript.

    bool noFragStartPosDist; // Don't learn a non-uniform start distribution

    bool noSeqBiasModel; // Don't learn and use a sequence-specific bias model.

    bool noRichEqClasses; // Don't use rich equivalence classes --- forget the
                          // aggregate weights for each transcript to each
                          // equivalence class of fragments.

    bool noFragLenFactor; // Don't account for the fragment length term in the
                          // likelihood.

    bool useReadCompat; // Give a fragment assignment a likelihood based on the compatibility
                        // between the manner in which it mapped and the expected read
                        // library format.

    double incompatPrior; // The prior probability that an alignment that disagrees with
                          // the provided library type could correspond to the true
                          // fragment origin.

    bool useErrorModel; // Learn and apply the error model when computing the likelihood
                        // of a given alignment.

    uint32_t numErrorBins; // Number of bins into which each read is divided
                           // when learning and applying the error model.

    double forgettingFactor; // The forgetting factor used to determine the
                             // learning schedule in the online inference algorithm.

    uint32_t numBurninFrags; // Number of mapped fragments required for burn-in

    uint32_t numPreBurninFrags; // Number of mapped fragments that are evaluated before
    				// applying the auxiliary models.

    uint32_t maxReadOccs; // Discard reads  mapping to more than this many places.

    uint32_t maxExpectedReadLen; // Maximum expected length of an observed read.

    bool useMassBanking; // DEPRECATED

    bool sensitive; // Perform splitting of long SMEMs into MEMs 
    
    bool extraSeedPass; // Perform extra pass trying to find seeds to cover the read

    bool disableMappingCache; // Don't write mapping results to temporary mapping cache file

    boost::filesystem::path outputDirectory; // Quant output directory

    boost::filesystem::path indexDirectory; // Index directory

    bool useVBOpt; // Use Variational Bayesian EM instead of "regular" EM in the batch passes
    bool useGSOpt; // Do Gibbs Sampling after optimization

    uint32_t numGibbsSamples; // Number of rounds of Gibbs sampling to perform 

    // Related to the fragment length distribution
    size_t fragLenDistMax;
    size_t fragLenDistPriorMean;
    size_t fragLenDistPriorSD;

    // Related to the logger
    std::shared_ptr<spdlog::logger> jointLog{nullptr};
    std::shared_ptr<spdlog::logger> fileLog{nullptr};

    // Related to caching and threading
    uint32_t mappingCacheMemoryLimit;
    uint32_t numThreads;
    uint32_t numQuantThreads;
    uint32_t numParseThreads;
};

#endif // SALMON_OPTS_HPP
