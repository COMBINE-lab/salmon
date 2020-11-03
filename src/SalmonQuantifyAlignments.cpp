
extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
}

// for cpp-format
#include "spdlog/fmt/fmt.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <tbb/concurrent_queue.h>

#include <boost/filesystem.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "AlignmentLibrary.hpp"
#include "BAMQueue.hpp"
#include "BAMUtils.hpp"
#include "ClusterForest.hpp"
#include "FASTAParser.hpp"
#include "LibraryFormat.hpp"
#include "MiniBatchInfo.hpp"
#include "SalmonExceptions.hpp"
#include "SalmonMath.hpp"
#include "Transcript.hpp"
#include "ProgramOptionsGenerator.hpp"

#include "SalmonDefaults.hpp"
#include "AlignmentModel.hpp"
#include "BiasParams.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "ForgettingMassCalculator.hpp"
#include "FragmentLengthDistribution.hpp"
#include "GZipWriter.hpp"
#include "NullFragmentFilter.hpp"
#include "ReadPair.hpp"
#include "SalmonConfig.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
#include "Sampler.hpp"
#include "TextBootstrapWriter.hpp"
#include "TranscriptCluster.hpp"
#include "spdlog/spdlog.h"
#include "pufferfish/Util.hpp"

namespace bfs = boost::filesystem;
using salmon::math::LOG_0;
using salmon::math::LOG_1;
using salmon::math::logAdd;
using salmon::math::logSub;

using MateStatus = pufferfish::util::MateStatus;

constexpr uint32_t miniBatchSize{1000};

template <typename FragT> using AlignmentBatch = std::vector<FragT>;

template <typename FragT>
using MiniBatchQueue = tbb::concurrent_queue<MiniBatchInfo<FragT>*>;

using PriorAbundanceVector = std::vector<double>;
using PosteriorAbundanceVector = std::vector<double>;

template <typename FragT>
using AlignmentLibraryT = AlignmentLibrary<FragT, EquivalenceClassBuilder<TGValue>>;

struct RefSeq {
  RefSeq(char* name, uint32_t len) : RefName(name), RefLength(len) {}
  std::string RefName;
  uint32_t RefLength;
};

/**
 * Multiple each element in the vector `vec` by the factor `scale`.
 */
template <typename T> void scaleBy(std::vector<T>& vec, T scale) {
  std::for_each(vec.begin(), vec.end(),
                [scale](T& ele) -> void { ele *= scale; });
}

/*
 * Tries _numTries_ times to get work from _workQueue_.  It returns
 * true immediately if it was able to find work, and false otherwise.
 */
template <typename FragT>
inline bool tryToGetWork(MiniBatchQueue<AlignmentGroup<FragT*>>& workQueue,
                         MiniBatchInfo<AlignmentGroup<FragT*>>*& miniBatch,
                         uint32_t numTries) {

  uint32_t attempts{1};
  bool foundWork = workQueue.try_pop(miniBatch);
  while (!foundWork and attempts < numTries) {
    foundWork = workQueue.try_pop(miniBatch);
    ++attempts;
  }
  return foundWork;
}


template <typename FragT>
void processMiniBatch(AlignmentLibraryT<FragT>& alnLib,
                      ForgettingMassCalculator& fmCalc,
                      uint64_t firstTimestepOfRound,
                      MiniBatchQueue<AlignmentGroup<FragT*>>& workQueue,
                      MiniBatchQueue<AlignmentGroup<FragT*>>* processedCache,
                      std::condition_variable& workAvailable,
                      std::mutex& cvmutex, volatile bool& doneParsing,
                      std::atomic<size_t>& activeBatches,
                      SalmonOpts& salmonOpts, BiasParams& observedBiasParams,
                      std::atomic<bool>& burnedIn, bool initialRound,
                      std::atomic<size_t>& processedReads) {

  // Seed with a real random value, if available
  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    std::random_device rd("/dev/urandom");
  #else
    std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 2020012

  auto& log = salmonOpts.jointLog;

  // Whether or not we are using "banking"
  bool useMassBanking = (!initialRound and salmonOpts.useMassBanking);
  double incompatPrior = salmonOpts.incompatPrior;
  bool useReadCompat = incompatPrior != salmon::math::LOG_1;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<> uni(
      0.0, 1.0 + std::numeric_limits<double>::min());

  // If we're auto detecting the library type
  auto* detector = alnLib.getDetector();
  bool autoDetect = (detector != nullptr) ? detector->isActive() : false;
  // If we haven't detected yet, nothing is incompatible
  if (autoDetect) {
    incompatPrior = salmon::math::LOG_1;
  }

  // EQClass
  auto& eqBuilder = alnLib.equivalenceClassBuilder();
  auto& readBiasFW = observedBiasParams.seqBiasModelFW;
  auto& readBiasRC = observedBiasParams.seqBiasModelRC;
  auto& observedGCMass = observedBiasParams.observedGCMass;
  auto& obsFwd = observedBiasParams.massFwd;
  auto& obsRC = observedBiasParams.massRC;
  auto& observedPosBiasFwd = observedBiasParams.posBiasFW;
  auto& observedPosBiasRC = observedBiasParams.posBiasRC;

  bool posBiasCorrect = salmonOpts.posBiasCorrect;
  bool gcBiasCorrect = salmonOpts.gcBiasCorrect;

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;
  using salmon::math::LOG_EPSILON;
  using salmon::math::logAdd;
  using salmon::math::logSub;

  auto orphanProb = salmonOpts.discardOrphansAln ? LOG_0 : LOG_EPSILON;

  // k-mers for sequence bias
  //Mer leftMer;
  //Mer rightMer;
  //Mer context;
  SBMer leftMer;
  SBMer rightMer;
  SBMer context;


  auto& refs = alnLib.transcripts();
  auto& clusterForest = alnLib.clusterForest();
  auto& fragmentQueue = alnLib.fragmentQueue();
  auto& alignmentGroupQueue = alnLib.alignmentGroupQueue();

  std::vector<FragmentStartPositionDistribution>& fragStartDists =
      alnLib.fragmentStartPositionDistributions();

  auto& fragLengthDist = *(alnLib.fragmentLengthDistribution());
  auto& alnMod = alnLib.alignmentModel();

  bool useFragLengthDist{!salmonOpts.noFragLengthDist};
  bool noFragLenFactor{salmonOpts.noFragLenFactor};

  double startingCumulativeMass =
      fmCalc.cumulativeLogMassAt(firstTimestepOfRound);
  LibraryFormat expectedLibraryFormat = alnLib.format();
  uint32_t numBurninFrags{salmonOpts.numBurninFrags};
  bool noLengthCorrection{salmonOpts.noLengthCorrection};

  uint32_t rangeFactorization{salmonOpts.rangeFactorizationBins};

  bool useAuxParams = (processedReads >= salmonOpts.numPreBurninFrags);
  bool singleEndLib = !alnLib.isPairedEnd();
  bool modelSingleFragProb = !salmonOpts.noSingleFragProb;
  size_t prevProcessedReads{0};
  size_t fragUpdateThresh{100000};

  distribution_utils::LogCMFCache logCMFCache(&fragLengthDist, singleEndLib);

  const size_t maxCacheLen{salmonOpts.fragLenDistMax};
  // Caches to avoid fld updates _within_ the set of alignments of a fragment 
  auto fetchPMF = [&fragLengthDist](size_t l) -> double { return fragLengthDist.pmf(l); };
  auto fetchCMF = [&fragLengthDist](size_t l) -> double { return fragLengthDist.cmf(l); };
  distribution_utils::IndexedVersionedCache<double> pmfCache(maxCacheLen);
  distribution_utils::IndexedVersionedCache<double> cmfCache(maxCacheLen);
  
  std::chrono::microseconds sleepTime(1);
  MiniBatchInfo<AlignmentGroup<FragT*>>* miniBatch = nullptr;
  bool updateCounts = initialRound;
  size_t numTranscripts = refs.size();

  double maxZeroFrac{0.0};

  auto isUnexpectedOrphan = [](FragT* aln, LibraryFormat expectedLibFormat) -> bool {
    return (expectedLibFormat.type == ReadType::PAIRED_END and
            !aln->isPaired());
  };

  auto alignerType = alnLib.getAlignerType();

  bool haveASTag{false};
  bool firstTagCheck{true};
  
  // If we are dealing with RapMap or Pufferfish mappings, then we 
  // don't expect CIGAR strings.  However, with recent versions 
  // at least any version that we could consider to be reasonable 
  // to run in alignment mode, we should have the `AS` tag. We set 
  // this the variable so we know this and deal with 
  // alignment score appropriately (and don't expect CIGAR strings).
  // NOTE: If we are reading in pufferfish or RapMap alignments, we expect the AS tag,
  // we will raise an exception if we don't have it.
  bool useASWithoutCIGAR = ((alignerType == salmon::bam_utils::AlignerDetails::PUFFERFISH) or 
                            (alignerType == salmon::bam_utils::AlignerDetails::RAPMAP));

  auto getAlignerAssignedScore =
    [&haveASTag, &firstTagCheck, alignerType](FragT* aln) -> double {
      if (firstTagCheck and !haveASTag) {
        char* tp = bam_aux_find(aln->getRead1(), "AS");
        haveASTag = (tp != NULL);
      }
      double score{LOG_0};
      if (haveASTag and (alignerType == salmon::bam_utils::AlignerDetails::BOWTIE2)) {
        uint8_t* tl = reinterpret_cast<uint8_t*>(bam_aux_find(aln->getRead1(), "AS"));
        auto locScore = (tl != NULL) ? bam_aux_i(tl) : LOG_1;
        if (aln->isPaired()) {
          uint8_t* tr = reinterpret_cast<uint8_t*>(bam_aux_find(aln->getRead2(), "AS"));
          locScore += (tr != NULL) ? bam_aux_i(tr) : LOG_1;
        }
        score = locScore;
      } else {
        score = LOG_1;
      }
      firstTagCheck = false;
      return score;
    };

  while (!doneParsing or !workQueue.empty()) {
    uint32_t zeroProbFrags{0};

    // Try up to numTries times to get work from the queue before
    // giving up and waiting on the condition variable
    constexpr uint32_t numTries = 100;
    bool foundWork = tryToGetWork(workQueue, miniBatch, numTries);

    // If work wasn't immediately available, then wait for it using
    // a condition variable to avoid burning CPU cycles for no reason.
    if (!foundWork) {
      std::unique_lock<std::mutex> l(cvmutex);
      workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() {
        return workQueue.try_pop(miniBatch) or doneParsing;
      });
    }

    uint64_t batchReads{0};

    // If we actually got some work
    if (miniBatch != nullptr) {
      expectedLibraryFormat = alnLib.format();
      useAuxParams = (processedReads >= salmonOpts.numPreBurninFrags);
      bool considerCondProb = (useAuxParams or burnedIn);
      ++activeBatches;
      batchReads = 0;
      zeroProbFrags = 0;

      // double logForgettingMass = fmCalc();
      double logForgettingMass{0.0};
      uint64_t currentMinibatchTimestep{0};
      // logForgettingMass and currentMinibatchTimestep are OUT parameters!
      fmCalc.getLogMassAndTimestep(logForgettingMass, currentMinibatchTimestep);
      miniBatch->logForgettingMass = logForgettingMass;

      std::vector<AlignmentGroup<FragT*>*>& alignmentGroups =
          *(miniBatch->alignments);

      using TranscriptID = size_t;
      using HitIDVector = std::vector<size_t>;
      using HitProbVector = std::vector<double>;

      // If we are going to attempt to model single mappings (part of a fragment)
      // then cache the FLD cumulative distribution for this mini-batch.  If
      // we are burned in or this is a single-end library, then the CMF is already
      // cached, so we don't need to worry about doing this work ourselves.
      if (modelSingleFragProb) {
        logCMFCache.refresh(processedReads.load(), burnedIn.load());
      }

      {
        bool checkedASTag = false;
        // Iterate over each group of alignments (a group consists of all
        // alignments reported for a single read).  Distribute the read's mass
        // proportionally dependent on the current
        for (auto& alnGroup : alignmentGroups) {
          pmfCache.increment_generation();
          cmfCache.increment_generation();

          // EQCLASS
          std::vector<uint32_t> txpIDs;
          std::vector<double> auxProbs;
          double auxDenom = salmon::math::LOG_0;

          // The alignments must be sorted by transcript id
          alnGroup->sortHits();

          double sumOfAlignProbs{LOG_0};

          // update the cluster-level properties
          bool transcriptUnique{true};
          auto firstTranscriptID =
              alnGroup->alignments().front()->transcriptID();
          std::unordered_set<size_t> observedTranscripts;

          int32_t bestAS = std::numeric_limits<int32_t>::min();

          // TODO: Use the AS tag to get the appropriate score for 
          // reads aligned with pufferfish or RapMap.
          
          // If we are using the AS tag without CIGAR strings
          if (useASWithoutCIGAR) {
            // Loop over all of the alignments to find the maximum score, which 
            // we will need to compute the alignment probability.
            for (auto& aln : alnGroup->alignments()) {
              if (!checkedASTag) {
                checkedASTag = true;
                if (! aln->haveASTag() ) {
                  log->critical("salmon is being run in alignment mode with a SAM/BAM file generated by {}, but the "
                  "alignment records do not seem to have the AS tag.  In this case, quantifying from alignments is not "
                  "recommended or supported.", salmon::bam_utils::to_string(alignerType));
                  spdlog::drop_all();
                  std::exit(1);
                }
                haveASTag = true;
              }

              int32_t score = aln->getAS();
              bestAS = (score > bestAS) ? score : bestAS;
            }
          }

          //double maxLogAlnScore{LOG_0};
          int sidx{0};
          for (auto& aln : alnGroup->alignments()) {
            auto transcriptID = aln->transcriptID();
            auto& transcript = refs[transcriptID];
            transcriptUnique =
                transcriptUnique and (transcriptID == firstTranscriptID);

            double refLength =
                transcript.RefLength > 0 ? transcript.RefLength : 1.0;
            auto flen = aln->fragLen();
            // If we have a properly-paired read then use the "pedantic"
            // definition here.
            if (aln->isPaired() and aln->isInward()) {
              flen = aln->fragLengthPedantic(refLength);
            }

            // The probability of drawing a fragment of this length;
            double logFragProb = LOG_1;
            // if we are modeling fragment probabilities for single-end mappings
            // and this is either a single-end library or an orphan.
            if (modelSingleFragProb and useFragLengthDist and (singleEndLib or isUnexpectedOrphan(aln, expectedLibraryFormat))) {
              logFragProb = logCMFCache.getAmbigFragLengthProb(aln->fwd(), aln->pos(), aln->readLen(), transcript.CompleteLength, burnedIn.load());
            } else if (isUnexpectedOrphan(aln, expectedLibraryFormat)) {
              // If we are expecting a paired-end library, and this is an orphan,
              // then logFragProb should be small
              logFragProb = orphanProb;
              if (logFragProb == LOG_0) {
                continue;
              }
            }

            if (flen > 0.0 and aln->isPaired() and useFragLengthDist and
                considerCondProb) {
              
              size_t fl = flen;
              double lenProb = pmfCache.get_or_update(fl, fetchPMF);

              if (burnedIn) {
                /* condition fragment length prob on txp length */
                size_t rlen = static_cast<size_t>(refLength);
                double refLengthCM = cmfCache.get_or_update(fl, fetchCMF);

                bool computeMass =
                    fl < refLength and !salmon::math::isLog0(refLengthCM);
                logFragProb = (computeMass) ? (lenProb - refLengthCM)
                                            : salmon::math::LOG_EPSILON;
                if (computeMass and refLengthCM < lenProb) {
                  // Threading is hard!  It's possible that an update to the PMF
                  // snuck in between when we asked to cache the CMF and when
                  // the "burnedIn" variable was last seen as false.
                  log->info("reference length = {}, CMF[refLen] = {}, fragLen "
                            "= {}, PMF[fragLen] = {}",
                            refLength, std::exp(refLengthCM), aln->fragLen(),
                            std::exp(lenProb));
                }
              } else if (useAuxParams) {
                logFragProb = lenProb;
              }
            }

            /*
            if (!salmonOpts.noFragLengthDist and useAuxParams) {
                if(aln->isPaired() and flen > 0) {
                    logFragProb = fragLengthDist.pmf(static_cast<size_t>(flen));
                }
            }
            */

            // TESTING
            if (noFragLenFactor) {
              logFragProb = LOG_1;
            }

            if (autoDetect) {
              detector->addSample(aln->libFormat());
              if (detector->canGuess()) {
                detector->mostLikelyType(alnLib.getFormat());
                expectedLibraryFormat = alnLib.format();
                incompatPrior = salmonOpts.incompatPrior;
                autoDetect = false;
              } else if (!detector->isActive()) {
                expectedLibraryFormat = alnLib.format();
                incompatPrior = salmonOpts.incompatPrior;
                autoDetect = false;
              }
            }
            // @TODO: handle this case better
            // double fragProb = cdf(fragLengthDist, fragLength + 0.5) -
            // cdf(fragLengthDist, fragLength - 0.5);  fragProb =
            // std::max(fragProb, 1e-3);  fragProb /= cdf(fragLengthDist,
            // refLength);

            // The alignment probability is the product of a
            // transcript-level term (based on abundance and) an
            // alignment-level term.
            double logRefLength{salmon::math::LOG_0};
            if (noLengthCorrection) {
              logRefLength = 1.0;
            } else if (salmonOpts.noEffectiveLengthCorrection or !burnedIn) {
              logRefLength = std::log(transcript.RefLength);
            } else {
              logRefLength = transcript.getCachedLogEffectiveLength();
            }

            // The probability that the fragments align to the given strands in
            // the given orientations.
            bool isCompat = salmon::utils::isCompatible(
                aln->libFormat(), expectedLibraryFormat, aln->pos(), aln->fwd(),
                aln->mateStatus());
            double logAlignCompatProb = isCompat ? LOG_1 : incompatPrior;
            if (!isCompat and salmonOpts.ignoreIncompat) {
              aln->logProb = salmon::math::LOG_0;
              continue;
            }
            /*
      double logAlignCompatProb =
          (useReadCompat) ?
          (salmon::utils::logAlignFormatProb(
                aln->libFormat(),
                expectedLibraryFormat,
                aln->pos(),
                aln->fwd(), aln->mateStatus(), salmonOpts.incompatPrior)
          ) : LOG_1;

      if (logAlignCompatProb != salmon::math::LOG_1) {
          aln->logProb = salmon::math::LOG_0;
          std::cerr <<"here\n";
          continue;
      }
            */

            // Adjustment to the likelihood due to the
            // error model
            double errLike = salmon::math::LOG_1;
            if (useASWithoutCIGAR) {
              int32_t alnScore = aln->getAS();
              // NOTE: we work directly in log space here, so the log
              // prob is -scoreExp * (S-w) rather than exp^(-scoreExp * (S-w))
              errLike = -salmonOpts.scoreExp * (bestAS - alnScore);
            } else if (useAuxParams and salmonOpts.useErrorModel) {
              errLike = alnMod.logLikelihood(*aln, transcript);
              ++sidx;
            }

            // Allow for a non-uniform fragment start position distribution
            double startPosProb{-logRefLength};
            if (aln->isPaired() and !noLengthCorrection) {
              startPosProb = (flen <= refLength)
                                 ? -std::log(refLength - flen + 1)
                                 : salmon::math::LOG_EPSILON;
            }

            double fragStartLogNumerator{salmon::math::LOG_1};
            double fragStartLogDenominator{salmon::math::LOG_1};

            auto hitPos = aln->left();

            // The total auxiliary probabilty is the product (sum in log-space)
            // of The fragment length probabilty The mapping score (under error
            // model) probability The fragment compatibility probability

            // The auxProb does *not* account for the start position
            // probability!
            double auxProb = logFragProb + errLike + logAlignCompatProb;

            // The overall mass of this transcript, which is used to
            // account for this transcript's relaive abundance
            double transcriptLogCount = transcript.mass(initialRound);

            if (transcriptLogCount != LOG_0 and auxProb != LOG_0 and
                startPosProb != LOG_0) {
              aln->logProb = transcriptLogCount + auxProb + startPosProb;

              sumOfAlignProbs = logAdd(sumOfAlignProbs, aln->logProb);
              if (updateCounts and observedTranscripts.find(transcriptID) ==
                                       observedTranscripts.end()) {
                refs[transcriptID].addTotalCount(1);
                observedTranscripts.insert(transcriptID);
              }
              // EQCLASS
              txpIDs.push_back(transcriptID);
              auxProbs.push_back(auxProb);
              auxDenom = salmon::math::logAdd(auxDenom, auxProb);

            } else {
              aln->logProb = LOG_0;
            }
          }

          // If we have a 0-probability fragment
          if (sumOfAlignProbs == LOG_0) {
            ++zeroProbFrags;
            ++batchReads;
            continue;
          }

          // EQCLASS
          double auxProbSum{0.0};
          for (auto& p : auxProbs) {
            p = std::exp(p - auxDenom);
            auxProbSum += p;
          }

          if (txpIDs.size() > 0) {

            if (rangeFactorization > 0) {
              int txpsSize = txpIDs.size();
              int rangeCount = std::sqrt(txpsSize) + rangeFactorization;

              for (int32_t i = 0; i < txpsSize; i++) {
                int rangeNumber = auxProbs[i] * rangeCount;
                txpIDs.push_back(rangeNumber);
              }
            }

            TranscriptGroup tg(txpIDs);
            eqBuilder.addGroup(std::move(tg), auxProbs);
          }

          // Are we doing bias correction?
          bool needBiasSample = salmonOpts.biasCorrect;

          // Normalize the scores
          for (auto& aln : alnGroup->alignments()) {
            if (aln->logProb == LOG_0) {
              continue;
            }
            aln->logProb -= sumOfAlignProbs;

            auto transcriptID = aln->transcriptID();
            auto& transcript = refs[transcriptID];

            double newMass = logForgettingMass + aln->logProb;
            transcript.addMass(newMass);

            // ---- Collect seq-specific bias samples ------ //
            auto getCIGARLength = [](bam_seq_t* s) -> uint32_t {
              auto cl = bam_cigar_len(s);
              decltype(cl) k, end;
              end = 0; // bam_pos(s);
              uint32_t* cigar = bam_cigar(s);
              for (k = 0; k < cl; ++k) {
                int op = cigar[k] & BAM_CIGAR_MASK;
                if (BAM_CONSUME_SEQ(op)) {
                  end += cigar[k] >> BAM_CIGAR_SHIFT;
                }
              }
              return end;
            };

            bool success = false;
            if (needBiasSample and salmonOpts.numBiasSamples > 0) {
              const char* txpStart = transcript.Sequence();
              const char* txpEnd = txpStart + transcript.RefLength;
              if (aln->isPaired()) {
                ReadPair* alnp = reinterpret_cast<ReadPair*>(aln);
                bam_seq_t* r1 = alnp->read1;
                bam_seq_t* r2 = alnp->read2;
                if (r1 != nullptr and r2 != nullptr) {
                  int32_t pos1 = bam_pos(r1);
                  bool fwd1{bam_strand(r1) == 0};
                  int32_t startPos1 =
                      fwd1 ? pos1 : pos1 + getCIGARLength(r1) - 1;

                  int32_t pos2 = bam_pos(r2);
                  bool fwd2{bam_strand(r2) == 0};
                  int32_t startPos2 =
                      fwd2 ? pos2 : pos2 + getCIGARLength(r2) - 1;

                  // Shouldn't be from the same strand and they should be in the
                  // right order
                  if ((fwd1 != fwd2) and // Shouldn't be from the same strand
                      (startPos1 > 0 and startPos1 < static_cast<int32_t>(transcript.RefLength)) and
                      (startPos2 > 0 and startPos2 < static_cast<int32_t>(transcript.RefLength))) {

                    const char* readStart1 = txpStart + startPos1;
                    auto& readBias1 = (fwd1) ? readBiasFW : readBiasRC;

                    const char* readStart2 = txpStart + startPos2;
                    auto& readBias2 = (fwd2) ? readBiasFW : readBiasRC;

                    int32_t fwPre = readBias1.contextBefore(!fwd1);
                    int32_t fwPost = readBias1.contextAfter(!fwd1);

                    int32_t rcPre = readBias2.contextBefore(!fwd2);
                    int32_t rcPost = readBias2.contextAfter(!fwd2);

                    bool read1RC = !fwd1;
                    bool read2RC = !fwd2;

                    if ((startPos1 >= readBias1.contextBefore(read1RC) and
                         startPos1 + readBias1.contextAfter(read1RC) <
                         static_cast<int32_t>(transcript.RefLength)) and
                        (startPos2 >= readBias2.contextBefore(read2RC) and
                         startPos2 + readBias2.contextAfter(read2RC) <
                         static_cast<int32_t>(transcript.RefLength))) {

                      int32_t fwPos = (fwd1) ? startPos1 : startPos2;
                      int32_t rcPos = (fwd1) ? startPos2 : startPos1;
                      if (fwPos < rcPos) {
                        leftMer.fromChars(txpStart + startPos1 -
                                           readBias1.contextBefore(read1RC));
                        rightMer.fromChars(txpStart + startPos2 -
                                            readBias2.contextBefore(read2RC));
 
                        if (read1RC) {
                          leftMer.rc();
                        } else {
                          rightMer.rc();
                        }

                        success = readBias1.addSequence(leftMer, 1.0);
                        success = readBias2.addSequence(rightMer, 1.0);
                      }
                    }
                  }
                }
              } else { // unpaired read
                UnpairedRead* alnp = reinterpret_cast<UnpairedRead*>(aln);
                bam_seq_t* r1 = alnp->read;
                if (r1 != nullptr) {
                  int32_t pos1 = bam_pos(r1);
                  bool fwd1{bam_strand(r1) == 0};
                  int32_t startPos1 =
                      fwd1 ? pos1 : pos1 + getCIGARLength(r1) - 1;

                  if (startPos1 > 0 and startPos1 < static_cast<int32_t>(transcript.RefLength)) {

                    const char* txpStart = transcript.Sequence();
                    const char* txpEnd = txpStart + transcript.RefLength;

                    const char* readStart1 = txpStart + startPos1;
                    auto& readBias1 = (fwd1) ? readBiasFW : readBiasRC;

                    if (startPos1 >= readBias1.contextBefore(!fwd1) and
                        startPos1 + readBias1.contextAfter(!fwd1) <
                        static_cast<int32_t>(transcript.RefLength)) {
                      context.fromChars(txpStart + startPos1 -
                                         readBias1.contextBefore(!fwd1));

                      if (!fwd1) {
                        context.rc();
                      }
                      success = readBias1.addSequence(context, 1.0);
                    }
                  }
                }
              } // end unpaired read
              if (success) {
                salmonOpts.numBiasSamples -= 1;
                needBiasSample = false;
              }
            }
            // ---- Collect seq-specific bias samples ------ //

            /**
             * Update the auxiliary models.
             **/
            // Paired-end
            if (aln->isPaired()) {
              // TODO: Is this right for *all* library types?
              if (aln->fwd()) {
                obsFwd = salmon::math::logAdd(obsFwd, aln->logProb);
              } else {
                obsRC = salmon::math::logAdd(obsRC, aln->logProb);
              }
            } else if (aln->libFormat().type == ReadType::SINGLE_END) {
              // Single-end or orphan
              if (aln->libFormat().strandedness == ReadStrandedness::S) {
                obsFwd = salmon::math::logAdd(obsFwd, aln->logProb);
              } else {
                obsRC = salmon::math::logAdd(obsRC, aln->logProb);
              }
            }

            if (posBiasCorrect) {
              auto lengthClassIndex = transcript.lengthClassIndex();
              switch (aln->mateStatus()) {
              case MateStatus::PAIRED_END_PAIRED: {
                // TODO: Handle the non opposite strand case
                if (aln->isInward()) {
                  auto* read1 = aln->getRead1();
                  auto* read2 = aln->getRead2();
                  int32_t posFW = aln->fwd() ? bam_pos(read1) : bam_pos(read2) + bam_seq_len(read2);
                  int32_t posRC = aln->fwd() ? bam_pos(read2) + bam_seq_len(read2) : bam_pos(read1);
                  posFW = posFW < 0 ? 0 : posFW;
                  posFW = posFW >= static_cast<int32_t>(transcript.RefLength) ?
                    static_cast<int32_t>(transcript.RefLength) - 1
                    : posFW;
                  posRC = posRC < 0 ? 0 : posRC;
                  posRC = posRC >= static_cast<int32_t>(transcript.RefLength) ?
                    static_cast<int32_t>(transcript.RefLength) - 1
                    : posRC;
                  observedPosBiasFwd[lengthClassIndex].addMass(
                                                               posFW, transcript.RefLength, aln->logProb);
                  observedPosBiasRC[lengthClassIndex].addMass(
                                                              posRC, transcript.RefLength, aln->logProb);
                }
              } break;
              case MateStatus::PAIRED_END_LEFT:
              case MateStatus::PAIRED_END_RIGHT:
              case MateStatus::SINGLE_END: {
                int32_t pos = aln->pos();
                pos = pos < 0 ? 0 : pos;
                pos = pos >= static_cast<int32_t>(transcript.RefLength) ?
                  static_cast<int32_t>(transcript.RefLength) - 1 : pos;
                if (aln->fwd()) {
                  observedPosBiasFwd[lengthClassIndex].addMass(
                                                               pos, transcript.RefLength, aln->logProb);
                } else {
                  observedPosBiasRC[lengthClassIndex].addMass(
                                                              pos, transcript.RefLength, aln->logProb);
                }
              } break;
              default:
                break;
              }
            }

            // Collect the GC-fragment bias samples
            if (gcBiasCorrect) {
              if (aln->isPaired()) {
                ReadPair* alnp = reinterpret_cast<ReadPair*>(aln);
                bam_seq_t* r1 = alnp->read1;
                bam_seq_t* r2 = alnp->read2;
                if (r1 != nullptr and r2 != nullptr) {
                  bool fwd1{bam_strand(r1) == 0};
                  bool fwd2{bam_strand(r2) == 0};
                  int32_t start = alnp->left();
                  int32_t stop = alnp->right();

                  if (start >= 0 and stop < static_cast<int32_t>(transcript.RefLength)) {
                    bool valid{false};
                    auto desc = transcript.gcDesc(start, stop, valid);
                    if (valid) {
                      observedGCMass.inc(desc, aln->logProb);
                    }
                  }
                }
              } else if (expectedLibraryFormat.type == ReadType::SINGLE_END) {
                // Both expected and observed should be single end here
                UnpairedRead* alnp = reinterpret_cast<UnpairedRead*>(aln);
                bam_seq_t* r = alnp->read;
                if (r != nullptr) {
                  bool fwd{alnp->fwd()};
                  // For single-end reads, simply assume that every fragment
                  // has a length equal to the conditional mean (given the
                  // current transcript's length).
                  auto cmeans = alnLib.condMeans();
                  auto cmean = static_cast<int32_t>(
                      (transcript.RefLength >= cmeans.size())
                          ? cmeans.back()
                          : cmeans[transcript.RefLength]);
                  int32_t start =
                      fwd ? alnp->pos() : std::max(0, alnp->pos() - cmean);
                  int32_t stop = start + cmean;
                  // WITH CONTEXT
                  if (start >= 0 and stop < static_cast<int32_t>(transcript.RefLength)) {
                    bool valid{false};
                    auto desc = transcript.gcDesc(start, stop, valid);
                    if (valid) {
                      observedGCMass.inc(desc, aln->logProb);
                    }
                  }
                }
              }
            }
            // END: GC-fragment bias

            double r = uni(eng);
            if (!burnedIn and r < std::exp(aln->logProb)) {
              /**
               * Update the bias sequence-specific bias model
               **/

              /*
              if (needBiasSample and salmonOpts.numBiasSamples > 0 and isPaired)
              {
              // the "start" position is the leftmost position if
              // we hit the forward strand, and the leftmost
              // position + the read length if we hit the reverse complement
              bam_seq_t* r = aln->get5PrimeRead();
              if (r) {
                  bool fwd{bam_strand(r) == 0};
                  int32_t pos{bam_pos(r)};
                  int32_t startPos = fwd ? pos : pos + bam_seq_len(r);
                  auto dir = salmon::utils::boolToDirection(fwd);

                                  if (startPos > 0 and startPos <
              transcript.RefLength) { auto& readBias = (fwd) ? readBiasFW :
              readBiasRC; const char* txpStart = transcript.Sequence(); const
              char* readStart = txpStart + startPos; const char* txpEnd =
              txpStart + transcript.RefLength; bool success =
              readBias.update(txpStart, readStart, txpEnd, dir); if (success) {
                                          salmonOpts.numBiasSamples -= 1;
                                          needBiasSample = false;
                                      }
                                  }
                              }
                          }
                          */

              // Update the error model
              if (!useASWithoutCIGAR and salmonOpts.useErrorModel) {
                auto alignerScore = getAlignerAssignedScore(aln);
                alnMod.update(*aln, transcript, alignerScore, logForgettingMass);
              }
              // Update the fragment length distribution
              if (aln->isPaired() and !salmonOpts.noFragLengthDist) {
                double fragLength =
                    aln->fragLengthPedantic(transcript.RefLength);
                if (fragLength > 0) {
                  fragLengthDist.addVal(fragLength, logForgettingMass);
                }
              }
            }
          }

          // update the single target transcript
          if (transcriptUnique) {
            if (updateCounts) {
              refs[firstTranscriptID].addUniqueCount(1);
            }
            clusterForest.updateCluster(firstTranscriptID, 1, logForgettingMass,
                                        updateCounts);
          } else { // or the appropriate clusters
            // ughh . . . C++ still has some very rough edges
            clusterForest.template mergeClusters<FragT>(
                alnGroup->alignments().begin(), alnGroup->alignments().end());
            clusterForest.updateCluster(
                alnGroup->alignments().front()->transcriptID(), 1,
                logForgettingMass, updateCounts);
          }

          ++batchReads;
        } // end read group
      }   // end timer

      double individualTotal = LOG_0;
      {
        /*
        // M-step
        for (auto kv = hitList.begin(); kv != hitList.end(); ++kv) {
            auto transcriptID = kv->first;
            // The target must be a valid transcript
            if (transcriptID >= numTranscripts or transcriptID < 0) {std::cerr
        << "index " << transcriptID << " out of bounds\n"; }

            auto& transcript = refs[transcriptID];

            // The prior probability
            double hitMass{LOG_0};

            // The set of alignments that match transcriptID
            auto& hits = kv->second;
            std::for_each(hits.begin(), hits.end(), [&](FragT* aln) -> void {
                    if (!std::isfinite(aln->logProb)) { log->warn("hitMass =
        {}\n", aln->logProb); } hitMass = logAdd(hitMass, aln->logProb);
            });

            // Lock the target
            if (hitMass == LOG_0) {
                log->warn("\n\n\n\nA set of *valid* alignments for a read
        appeared to " "have 0 probability.  This should not happen.  Please
        report " "this bug.  exiting!");

                std::cerr << "\n\n\n\nA set of *valid* alignments for a read
        appeared to "
                          << "have 0 probability.  This should not happen.
        Please report "
                          << "this bug.  exiting!";
                std::exit(1);
            }

            double updateMass = logForgettingMass + hitMass;
            individualTotal = logAdd(individualTotal, updateMass);
            transcript.addMass(updateMass);

           // unlock the target
        } // end for
        */
      } // end timer

      // If we're not keeping around a cache, then
      // reclaim the memory for these fragments and alignments
      // and delete the mini batch.
      if (processedCache == nullptr) {
        miniBatch->release(fragmentQueue, alignmentGroupQueue);
        delete miniBatch;
      } else {
        // Otherwise, just put the mini-batch on the processed queue
        // to be re-used in the next round.
        processedCache->push(miniBatch);
      }
      --activeBatches;
      processedReads += batchReads;
      if (processedReads >= numBurninFrags and !burnedIn) {
        // NOTE: only one thread should succeed here, and that
        // thread will set burnedIn to true
        alnLib.updateTranscriptLengthsAtomic(burnedIn);
        fragLengthDist.cacheCMF();
      }

      if (zeroProbFrags > 0) {
        maxZeroFrac =
            std::max(maxZeroFrac,
                     static_cast<double>(100.0 * zeroProbFrags) / batchReads);
      }
    }

    miniBatch = nullptr;
  } // nothing left to process

  if (maxZeroFrac > 0.0) {
    log->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero "
              "probability fragments",
              maxZeroFrac);
  }
}

/**
 *  Quantify the targets given in the file `transcriptFile` using the
 *  alignments given in the file `alignmentFile`, and write the results
 *  to the file `outputFile`.  The reads are assumed to be in the format
 *  specified by `libFmt`.
 *
 */
template <typename FragT>
bool quantifyLibrary(AlignmentLibraryT<FragT>& alnLib,
                     size_t numRequiredFragments, SalmonOpts& salmonOpts) {

  std::atomic<bool> burnedIn{salmonOpts.numBurninFrags == 0};

  auto& refs = alnLib.transcripts();
  size_t numTranscripts = refs.size();
  size_t numObservedFragments{0};

  auto& fileLog = salmonOpts.fileLog;
  bool useMassBanking = salmonOpts.useMassBanking;

  MiniBatchQueue<AlignmentGroup<FragT*>> workQueue;
  MiniBatchQueue<AlignmentGroup<FragT*>>* workQueuePtr{&workQueue};
  MiniBatchQueue<AlignmentGroup<FragT*>> processedCache;
  MiniBatchQueue<AlignmentGroup<FragT*>>* processedCachePtr{nullptr};

  ForgettingMassCalculator fmCalc(salmonOpts.forgettingFactor);

  // Up-front, prefill the forgetting mass schedule for up to
  // 1 billion reads
  size_t prefillSize = (1000000000) / miniBatchSize;
  fmCalc.prefill(prefillSize);

  size_t batchNum{0};
  std::atomic<size_t> totalProcessedReads{0};
  bool initialRound{true};
  bool haveCache{false};
  bool doReset{true};
  bool posBiasCorrect{salmonOpts.posBiasCorrect};
  bool gcBiasCorrect{salmonOpts.gcBiasCorrect};
  size_t maxCacheSize{salmonOpts.mappingCacheMemoryLimit};

  NullFragmentFilter<FragT>* nff = nullptr;
  bool terminate{false};

  // Give ourselves some space
  fmt::print(stderr, "\n\n\n\n");

  while (numObservedFragments < numRequiredFragments and !terminate) {
    if (!initialRound) {

      size_t numToCache = (useMassBanking)
                              ? (alnLib.numMappedFragments() -
                                 alnLib.numUniquelyMappedFragments())
                              : (alnLib.numMappedFragments());

      if (haveCache) {
        std::swap(workQueuePtr, processedCachePtr);
        doReset = false;
        fmt::print(stderr, "\n\n");
      } else if (numToCache <= maxCacheSize) {
        processedCachePtr = &processedCache;
        doReset = true;
        fmt::print(stderr, "\n");
      }

      if (doReset and
          !alnLib.reset(
              true,          /* increment # of passes */
              nff,           /* fragment filter */
              useMassBanking /* only process ambiguously mapped fragments*/
              )) {
        fmt::print(
            stderr,
            "\n\n======== WARNING ========\n"
            "A provided alignment file "
            "is not a regular file and therefore can't be read from "
            "more than once.\n\n"
            "We observed only {} fragments when we wanted at least {}.\n\n"
            "Please consider re-running Salmon with these alignments "
            "as a regular file!\n"
            "==========================\n\n",
            numObservedFragments, numRequiredFragments);
        break;
      }
    }

    volatile bool doneParsing{false};
    std::condition_variable workAvailable;
    std::mutex cvmutex;
    std::vector<std::thread> workers;
    std::atomic<size_t> activeBatches{0};
    auto currentQuantThreads =
        (haveCache) ? salmonOpts.numQuantThreads + salmonOpts.numParseThreads
                    : salmonOpts.numQuantThreads;

    uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
    if (firstTimestepOfRound > 1) {
      firstTimestepOfRound -= 1;
    }

    /** sequence-specific and GC-fragment bias vectors --- each thread gets it's
     * own **/
    std::vector<BiasParams> observedBiasParams(
        currentQuantThreads, BiasParams(salmonOpts.numConditionalGCBins,
                                        salmonOpts.numFragGCBins, false));

    for (uint32_t i = 0; i < currentQuantThreads; ++i) {
      workers.emplace_back(
          processMiniBatch<FragT>, std::ref(alnLib), std::ref(fmCalc),
          firstTimestepOfRound, std::ref(*workQueuePtr), processedCachePtr,
          std::ref(workAvailable), std::ref(cvmutex), std::ref(doneParsing),
          std::ref(activeBatches), std::ref(salmonOpts),
          std::ref(observedBiasParams[i]), std::ref(burnedIn), initialRound,
          std::ref(totalProcessedReads));
    }

    if (!haveCache) {
      size_t numProc{0};

      BAMQueue<FragT>& bq = alnLib.getAlignmentGroupQueue();
      std::vector<AlignmentGroup<FragT*>*>* alignments =
          new std::vector<AlignmentGroup<FragT*>*>;
      alignments->reserve(miniBatchSize);
      AlignmentGroup<FragT*>* ag;

      bool alignmentGroupsRemain = bq.getAlignmentGroup(ag);
      while (alignmentGroupsRemain or alignments->size() > 0) {
        if (alignmentGroupsRemain) {
          alignments->push_back(ag);
        }
        // If this minibatch has reached the size limit, or we have nothing
        // left to fill it up with
        if (alignments->size() >= miniBatchSize or !alignmentGroupsRemain) {
          ++batchNum;
          double logForgettingMass = 0.0;
          MiniBatchInfo<AlignmentGroup<FragT*>>* mbi =
              new MiniBatchInfo<AlignmentGroup<FragT*>>(batchNum, alignments,
                                                        logForgettingMass);
          workQueuePtr->push(mbi);
          {
            std::unique_lock<std::mutex> l(cvmutex);
            workAvailable.notify_one();
          }
          alignments = new std::vector<AlignmentGroup<FragT*>*>;
          alignments->reserve(miniBatchSize);
        }

        if ((numProc % 1000000 == 0) or !alignmentGroupsRemain) {
          fmt::print(stderr, "\r\r{}processed{} {} {}reads in current round{}",
                     ioutils::SET_GREEN, ioutils::SET_RED, numProc,
                     ioutils::SET_GREEN, ioutils::RESET_COLOR);
          fileLog->info("quantification processed {} fragments so far\n",
                        numProc);
        }

        ++numProc;
        alignmentGroupsRemain = bq.getAlignmentGroup(ag);
      }
      fmt::print(stderr, "\n");

      // Free the alignments and the vector holding them
      if (processedCachePtr == nullptr) {
        for (auto& alnGroup : *alignments) {
          alnGroup->alignments().clear();
          delete alnGroup;
          alnGroup = nullptr;
        }
        delete alignments;
      }
    } else {
      fmt::print(stderr, "\n");
    }

    doneParsing = true;

    /**
     * This could be a problem for small sets of alignments --- make sure the
     * work queue is empty!!
     * --- Thanks for finding a dataset that exposes this bug, Richard
     * (Smith-Unna)!
     */
    size_t t = 0;
    while (!workQueuePtr->empty()) {
      std::unique_lock<std::mutex> l(cvmutex);
      workAvailable.notify_one();
    }

    size_t tnum{0};
    for (auto& t : workers) {
      fmt::print(stderr, "\r\rkilling thread {} . . . ", tnum++);
      {
        std::unique_lock<std::mutex> l(cvmutex);
        workAvailable.notify_all();
      }
      t.join();
      fmt::print(stderr, "done");
    }
    fmt::print(stderr, "\n\n");

    numObservedFragments += alnLib.numMappedFragments();

    // If we don't have a sufficient number of mapped fragments, then
    // complain here!

    // If we don't have a sufficient number of assigned fragments, then
    // complain here!
    if (numObservedFragments < salmonOpts.minRequiredFrags) {
      throw InsufficientAssignedFragments(numObservedFragments,
                                          salmonOpts.minRequiredFrags);
    }

    /**
     *
     * Aggregate thread-local bias parameters
     *
     **/
    // Set the global distribution based on the sum of local
    // distributions.
    double gcFracFwd{0.0};
    double globalMass{salmon::math::LOG_0};
    double globalFwdMass{salmon::math::LOG_0};
    auto& globalGCMass = alnLib.observedGC();
    for (auto& gcp : observedBiasParams) {
      auto& gcm = gcp.observedGCMass;
      globalGCMass.combineCounts(gcm);

      auto& fw =
          alnLib.readBiasModelObserved(salmon::utils::Direction::FORWARD);
      auto& rc = alnLib.readBiasModelObserved(
          salmon::utils::Direction::REVERSE_COMPLEMENT);

      auto& fwloc = gcp.seqBiasModelFW;
      auto& rcloc = gcp.seqBiasModelRC;
      fw.combineCounts(fwloc);
      rc.combineCounts(rcloc);

      /**
       * positional biases
       **/
      auto& posBiasesFW = alnLib.posBias(salmon::utils::Direction::FORWARD);
      auto& posBiasesRC =
        alnLib.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      for (size_t i = 0; i < posBiasesFW.size(); ++i) {
        posBiasesFW[i].combine(gcp.posBiasFW[i]);
        posBiasesRC[i].combine(gcp.posBiasRC[i]);
      }

      globalMass = salmon::math::logAdd(globalMass, gcp.massFwd);
      globalMass = salmon::math::logAdd(globalMass, gcp.massRC);
      globalFwdMass = salmon::math::logAdd(globalFwdMass, gcp.massFwd);
    }
    globalGCMass.normalize();

    if (globalMass != salmon::math::LOG_0) {
      if (globalFwdMass != salmon::math::LOG_0) {
        gcFracFwd = std::exp(globalFwdMass - globalMass);
      }
      alnLib.setGCFracForward(gcFracFwd);
    }

    // finalize the positional biases
    if (salmonOpts.posBiasCorrect) {
      auto& posBiasesFW = alnLib.posBias(salmon::utils::Direction::FORWARD);
      auto& posBiasesRC =
        alnLib.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      for (size_t i = 0; i < posBiasesFW.size(); ++i) {
        posBiasesFW[i].finalize();
        posBiasesRC[i].finalize();
      }
    }
    /** END: aggregate thread-local bias parameters **/

    fmt::print(
        stderr,
        "# observed = {} / # required = {}\033[A\033[A\033[A\033[A\033[A",
        numObservedFragments, numRequiredFragments);

    if (initialRound) {
      salmonOpts.jointLog->info(
          "\n\n\nCompleted first pass through the alignment file.\n"
          "Total # of mapped reads : {}\n"
          "# of uniquely mapped reads : {}\n"
          "# ambiguously mapped reads : {}\n\n\n",
          alnLib.numMappedFragments(), alnLib.numUniquelyMappedFragments(),
          alnLib.numMappedFragments() - alnLib.numUniquelyMappedFragments());
    }

    initialRound = false;

    // If we're done our second pass and we've decided to use
    // the in-memory cache, then activate it now.
    if (!initialRound and processedCachePtr != nullptr) {
      haveCache = true;
    }
    // EQCLASS
    bool done = alnLib.equivalenceClassBuilder().finish();
    // skip the extra online rounds
    terminate = true;
    // END EQCLASS
  }

  fmt::print(stderr, "\n\n\n\n");

  // If we didn't achieve burnin, then at least compute effective
  // lengths and mention this to the user.
  if (alnLib.numMappedFragments() < salmonOpts.numBurninFrags) {
    std::atomic<bool> dummyBool{false};
    alnLib.updateTranscriptLengthsAtomic(dummyBool);
    salmonOpts.jointLog->warn("Only {} fragments were mapped, but the number "
                              "of burn-in fragments was set to {}.\n"
                              "The effective lengths have been computed using "
                              "the observed mappings.\n",
                              alnLib.numMappedFragments(),
                              salmonOpts.numBurninFrags);

  }

  // In this case, we have to give the structures held
  // in the cache back to the appropriate queues
  if (haveCache) {
    auto& fragmentQueue = alnLib.fragmentQueue();
    auto& alignmentGroupQueue = alnLib.alignmentGroupQueue();

    MiniBatchInfo<AlignmentGroup<FragT*>>* mbi = nullptr;
    while (!processedCache.empty()) {
      while (processedCache.try_pop(mbi)) {
        mbi->release(fragmentQueue, alignmentGroupQueue);
        delete mbi;
      }
    }
  }

  return burnedIn.load();
}

template <typename ReadT>
bool processSample(AlignmentLibraryT<ReadT>& alnLib, size_t requiredObservations,
                   SalmonOpts& sopt, boost::filesystem::path outputDirectory) {

  auto& jointLog = sopt.jointLog;
  // EQCLASS
  alnLib.equivalenceClassBuilder().setMaxResizeThreads(sopt.maxHashResizeThreads);
  alnLib.equivalenceClassBuilder().start();

  bool burnedIn = false;
  try {
    // if this is a single-end library, then the fld won't change
    if (!alnLib.isPairedEnd()) {
      alnLib.fragmentLengthDistribution()->cacheCMF();
    }
    burnedIn = quantifyLibrary<ReadT>(alnLib, requiredObservations, sopt);
  } catch (const InsufficientAssignedFragments& iaf) {
    jointLog->warn(iaf.what());
    GZipWriter gzw(outputDirectory, jointLog);
    gzw.writeEmptyAbundances(sopt, alnLib);
    // Write meta-information about the run
    std::vector<std::string> errors{"insufficient_assigned_fragments"};
    sopt.runStopTime = salmon::utils::getCurrentTimeAsString();
    gzw.writeEmptyMeta(sopt, alnLib, errors);
    std::exit(1);
  }

  GZipWriter gzw(outputDirectory, jointLog);

  if (!sopt.skipQuant) {
    // NOTE: A side-effect of calling the optimizer is that
    // the `EffectiveLength` field of each transcript is
    // set to its final value.
    CollapsedEMOptimizer optimizer;
    jointLog->info("starting optimizer");
    salmon::utils::normalizeAlphas(sopt, alnLib);
    bool optSuccess = optimizer.optimize(alnLib, sopt, 0.01, 10000);
    // If the optimizer didn't work, then bail out here.
    if (!optSuccess) {
      return false;
    }
    jointLog->info("finished optimizer");

    jointLog->info("writing output");
    // Write the main results
    gzw.writeAbundances(sopt, alnLib, false);

    if (sopt.numGibbsSamples > 0) {

      jointLog->info("Starting Gibbs Sampler");
      CollapsedGibbsSampler sampler;
      gzw.setSamplingPath(sopt);
      // The function we'll use as a callback to write samples
      std::function<bool(const std::vector<double>&)> bsWriter =
        [&gzw](const std::vector<double>& alphas) -> bool {
          return gzw.writeBootstrap(alphas, true);
        };

      bool sampleSuccess =
        sampler.sample(alnLib, sopt, bsWriter, sopt.numGibbsSamples);
      if (!sampleSuccess) {
        jointLog->error("Encountered error during Gibb sampling .\n"
                        "This should not happen.\n"
                        "Please file a bug report on GitHub.\n");
        return false;
      }
      jointLog->info("Finished Gibbs Sampler");
    } else if (sopt.numBootstraps > 0) {
      // The function we'll use as a callback to write samples
      std::function<bool(const std::vector<double>&)> bsWriter =
        [&gzw](const std::vector<double>& alphas) -> bool {
          return gzw.writeBootstrap(alphas);
        };

      jointLog->info("Staring Bootstrapping");
      gzw.setSamplingPath(sopt);
      bool bootstrapSuccess =
        optimizer.gatherBootstraps(alnLib, sopt, bsWriter, 0.01, 10000);
      jointLog->info("Finished Bootstrapping");
      if (!bootstrapSuccess) {
        jointLog->error("Encountered error during bootstrapping.\n"
                        "This should not happen.\n"
                        "Please file a bug report on GitHub.\n");
        return false;
      }
    }

    // bfs::path libCountFilePath = outputDirectory / "lib_format_counts.json";
    // alnLib.summarizeLibraryTypeCounts(libCountFilePath);

    if (sopt.sampleOutput) {
      // In this case, we should "re-convert" transcript
      // masses to be counts in log space
      auto nr = alnLib.numMappedFragments();
      for (auto& t : alnLib.transcripts()) {
        double m = t.mass(false) * nr;
        if (m > 0.0) {
          t.setMass(std::log(m));
        }
      }

      bfs::path sampleFilePath = outputDirectory / "postSample.bam";
      bool didSample = salmon::sampler::sampleLibrary<ReadT>(
                                                             alnLib, sopt, burnedIn, sampleFilePath, sopt.sampleUnaligned);
      if (!didSample) {
        jointLog->warn("There may have been a problem generating the sampled "
                       "output file; please check the log\n");
      }
    }
  } else if (sopt.dumpEqWeights) { // sopt.skipQuant == true
    jointLog->info("Finalizing combined weights for equivalence classes.");
      // if we are skipping the quantification, and we are dumping equivalence class weights,
      // then fill in the combinedWeights of the equivalence classes so that `--dumpEqWeights` makes sense.
      auto& eqVec =
        alnLib.equivalenceClassBuilder().eqVec();
      bool noRichEq = sopt.noRichEqClasses;
      bool useEffectiveLengths = !sopt.noEffectiveLengthCorrection;
      std::vector<Transcript>& transcripts = alnLib.transcripts();
      Eigen::VectorXd effLens(transcripts.size());

      for (size_t i = 0; i < transcripts.size(); ++i) {
        auto& txp = transcripts[i];
        effLens(i) = useEffectiveLengths
          ? std::exp(txp.getCachedLogEffectiveLength())
          : txp.RefLength;
      }

      for (size_t eqID = 0; eqID < eqVec.size(); ++eqID){
        // The vector entry
        auto& kv = eqVec[eqID];
        // The label of the equivalence class
        const TranscriptGroup& k = kv.first;
        // The size of the label
        size_t classSize = kv.second.weights.size(); // k.txps.size();
        // The weights of the label
        auto& v = kv.second;

        // Iterate over each weight and set it
        double wsum{0.0};

        for (size_t i = 0; i < classSize; ++i) {
          auto tid = k.txps[i];
          double el = effLens(tid);
          if (el <= 1.0) {
            el = 1.0;
          }
          if (noRichEq) {
            // Keep length factor separate for the time being
            v.weights[i] = 1.0;
          }
          // meaningful values.
          auto probStartPos = 1.0 / el;

          // combined weight
          double wt = sopt.eqClassMode ? v.weights[i] : v.count * v.weights[i] * probStartPos;
          v.combinedWeights.push_back(wt);
          wsum += wt;
        }

        double wnorm = 1.0 / wsum;
        for (size_t i = 0; i < classSize; ++i) {
          v.combinedWeights[i] = v.combinedWeights[i] * wnorm;
        }
      }
      jointLog->info("done.");
    }


  // If we are dumping the equivalence classes, then
  // do it here.
  if (sopt.dumpEq) {
    gzw.writeEquivCounts(sopt, alnLib);
  }

  sopt.runStopTime = salmon::utils::getCurrentTimeAsString();
  // Write meta-information about the run
  MappingStatistics mstats;
  gzw.writeMeta(sopt, alnLib, mstats);

  return true;
}

bool processEqClasses( AlignmentLibraryT<UnpairedRead>& alnLib, SalmonOpts& sopt,
                       boost::filesystem::path outputDirectory) {
  auto& jointLog = sopt.jointLog;
  GZipWriter gzw(outputDirectory, jointLog);

  // NOTE: A side-effect of calling the optimizer is that
  // the `EffectiveLength` field of each transcript is
  // set to its final value.
  CollapsedEMOptimizer optimizer;
  jointLog->info("starting optimizer");

  {
    // setting no effective length correction, as we are taking effective lens as input
    sopt.initUniform = true;
    sopt.eqClassMode = true;
    jointLog->warn("Using Uniform Prior");
  }

  //salmon::utils::normalizeAlphas(sopt, alnLib);
  bool optSuccess = optimizer.optimize(alnLib, sopt, 0.01, 10000);
  // If the optimizer didn't work, then bail out here.
  if (!optSuccess) {
    return false;
  }
  jointLog->info("finished optimizer");
  jointLog->info("writing output");
  jointLog->flush();

  // Write the main results
  gzw.writeAbundances(sopt, alnLib, true);
  sopt.runStopTime = salmon::utils::getCurrentTimeAsString();

  return true;
}

int salmonAlignmentQuantify(int argc, const char* argv[]) {
  using std::cerr;
  using std::vector;
  using std::string;
  namespace po = boost::program_options;
  namespace bfs = boost::filesystem;

  SalmonOpts sopt;
  sopt.numThreads = salmon::defaults::numThreads;

  size_t requiredObservations{50000000};
  int32_t numBiasSamples{0};
  salmon::ProgramOptionsGenerator pogen;

  auto inputOpt = pogen.getAlignmentInputOptions(sopt);
  auto basicOpt = pogen.getBasicOptions(sopt);
  auto alnSpecOpt = pogen.getAlignmentSpecificOptions(sopt);
  auto advancedOpt = pogen.getAdvancedOptions(numBiasSamples, sopt);
  auto hiddenOpt = pogen.getHiddenOptions(sopt);
  auto testingOpt = pogen.getTestingOptions(sopt);
  auto deprecatedOpt = pogen.getDeprecatedOptions(sopt);

  po::options_description all("salmon quant options");
  all.add(inputOpt).add(basicOpt).add(alnSpecOpt).add(advancedOpt).add(testingOpt).add(hiddenOpt).add(deprecatedOpt);

  po::options_description visible("salmon quant options");
  visible.add(inputOpt).add(basicOpt).add(alnSpecOpt).add(advancedOpt);

  po::variables_map vm;
  try {
    auto orderedOptions =
        po::command_line_parser(argc, argv).options(all).run();

    po::store(orderedOptions, vm);

    if (vm.count("help")) {
      auto hstring = R"(
Quant
==========
Perform dual-phase, alignment-based estimation of
transcript abundance from RNA-seq reads
)";
      std::cout << hstring << std::endl;
      std::cout << visible << std::endl;
      std::exit(0);
    }

    po::notify(vm);

    if (sopt.numThreads < 2) {
      fmt::print(stderr, "salmon requires at least 2 threads --- "
                         "setting # of threads = 2\n");
      sopt.numThreads = 2;
    }
    auto numThreads = sopt.numThreads;

    bool hasAlignments {false};
    bool hasEqclasses {false};
    vector<string> alignmentFileNames;
    string eqclassesFileName;
    if (vm.count("alignments")) {
      alignmentFileNames = vm["alignments"].as<vector<string>>();
      hasAlignments = true;
    }
    if (vm.count("eqclasses")) {
      eqclassesFileName = vm["eqclasses"].as<string>();
      hasEqclasses = true;
    }

    if ( !hasAlignments and !hasEqclasses ) {
      fmt::print(stderr, "salmon requires at least one alignment input\n"
                 "Neither alignments (BAM) nor eqclasses given as input. \n");
      std::exit(1);
    }

    if (sopt.forgettingFactor <= 0.5 or sopt.forgettingFactor > 1.0) {
      fmt::print(stderr,
                 "The forgetting factor must be in (0.5, 1.0], "
                 "but the value {} was provided\n",
                 sopt.forgettingFactor);
      std::exit(1);
    }

    std::stringstream commentStream;
    commentStream << "# salmon (alignment-based) v" << salmon::version << "\n";
    commentStream << "# [ program ] => salmon \n";
    commentStream << "# [ command ] => quant \n";
    for (auto& opt : orderedOptions.options) {
      commentStream << "# [ " << opt.string_key << " ] => {";
      for (auto& val : opt.value) {
        commentStream << " " << val;
      }
      commentStream << " }\n";
    }
    std::string commentString = commentStream.str();
    if (!sopt.quiet) {
      fmt::print(stderr, "{}", commentString);
    }

    sopt.alnMode = true;
    sopt.quantMode = SalmonQuantMode::ALIGN;
    bool optionsOK =
      salmon::utils::processQuantOptions(sopt, vm, numBiasSamples);
    if (!optionsOK) {
      if (sopt.jointLog) {
        sopt.jointLog->flush();
        spdlog::drop_all();
      }
      std::exit(1);
    }

    auto fileLog = sopt.fileLog;
    auto jointLog = sopt.jointLog;
    auto indexDirectory = sopt.indexDirectory;
    auto outputDirectory = sopt.outputDirectory;

    // ==== Library format processing ===
    vector<bfs::path> alignmentFiles;

    if ( hasAlignments ) {
      for (auto& alignmentFileName : alignmentFileNames) {
        bfs::path alignmentFile(alignmentFileName);
        if (!bfs::exists(alignmentFile)) {
          std::stringstream ss;
          ss << "The provided alignment file: " << alignmentFile
             << " does not exist!\n";
          throw std::invalid_argument(ss.str());
        } else {
          alignmentFiles.push_back(alignmentFile);
        }
      }
    } else {
      bfs::path alignmentFile(eqclassesFileName);
      if (!bfs::exists(alignmentFile)) {
        std::stringstream ss;
        ss << "The provided eqclasses file: " << alignmentFile
           << " does not exist!\n";
        throw std::invalid_argument(ss.str());
      } else {
        alignmentFiles.push_back(alignmentFile);
      }
    }

    // Just so we have the variable around
    LibraryFormat libFmt(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                         ReadStrandedness::U);
    // Get the library format string
    std::string libFmtStr = vm["libType"].as<std::string>();
    bool autoDetectFmt =
      (libFmtStr == "a" or
       libFmtStr == "A"); //(autoTypes.find(libFmtStr) != autoTypes.end());

    if ( hasEqclasses ) {
      libFmt = LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                             ReadStrandedness::U);
    } else {
      // If we're auto-detecting, set things up appropriately
      if (autoDetectFmt) {

        bool isPairedEnd = salmon::utils::peekBAMIsPaired(alignmentFiles.front());

        if (isPairedEnd) {
          libFmt = LibraryFormat(ReadType::PAIRED_END, ReadOrientation::TOWARD,
                                 ReadStrandedness::U);
        } else {
          libFmt = LibraryFormat(ReadType::SINGLE_END, ReadOrientation::NONE,
                                 ReadStrandedness::U);
        }
      } else { // Parse the provided type
        libFmt = salmon::utils::parseLibraryFormatStringNew(libFmtStr);
      }
    }

    if (libFmt.check()) {
      std::cerr << libFmt << "\n";
    } else {
      std::stringstream ss;
      ss << libFmt << " is invalid!";
      throw std::invalid_argument(ss.str());
    }
    // ==== END: Library format processing ===

    // The transcript file contains the target sequences
    bfs::path transcriptFile;
    if (hasAlignments) { transcriptFile = vm["targets"].as<std::string>(); }

    // Currently, one thread is used for parsing the alignment file.
    // Hopefully, in the future, samtools will implemented multi-threaded
    // BAM/SAM parsing, as this is the current bottleneck.  For the time
    // being, however, the number of quantification threads is the
    // total number of threads - 1.
    uint32_t numParseThreads =
        std::min(uint32_t(6),
                 std::max(uint32_t(2), uint32_t(std::ceil(numThreads / 2.0))));
    numThreads = std::max(numThreads, numParseThreads);
    uint32_t numQuantThreads =
        std::max(uint32_t(2), uint32_t(numThreads - numParseThreads));
    sopt.numQuantThreads = numQuantThreads;
    sopt.numParseThreads = numParseThreads;
    jointLog->info("numQuantThreads = {}", numQuantThreads);

    // Write out information about the command / run
    salmon::utils::writeCmdInfo(sopt, orderedOptions);

    bool success{false};

    switch (libFmt.type) {
    case ReadType::SINGLE_END: {
      // We can only do fragment GC bias correction, for the time being, with
      // paired-end reads
      if (sopt.gcBiasCorrect) {
        jointLog->warn(
            "Fragment GC bias correction is currently *experimental* "
            "in single-end libraries.  Please use this option "
            "with caution.");
        // sopt.gcBiasCorrect = false;
      }

      if ( hasEqclasses ) {
        std::vector<string> tnames;
        std::vector<double> tefflens;
        std::vector<uint32_t> eqclass_counts;
        std::vector<std::vector<uint32_t>> eqclasses;
        std::vector<std::vector<double>> auxs_vals;
        {
          // reading eqclass
          bool parseOK = salmon::utils::readEquivCounts(alignmentFiles[0], tnames, tefflens,
                                                        eqclasses, auxs_vals, eqclass_counts);
          if (!parseOK){
            jointLog->error("Eqclass Parsing error");
            exit(1);
          }

          std::stringstream errfmt;
          errfmt << "Found total " << eqclasses.size() << " eqclasses and "
                 << tnames.size() << " transcripts";
          jointLog->info(errfmt.str());
          jointLog->flush();

          if ( tefflens.size() == 0 ) {
            tefflens.resize(tnames.size(), 100);
            jointLog->warn("No effective lens found in the eqclass file;"
                           "Ignore this warning if using uniform prior");
          }

          if ( /*tefflens.size() != 0 and*/ (tefflens.size() != tnames.size()) ){
            std::stringstream errfmt;
            errfmt << "Number of effective lens: " << tefflens.size()
                   << " is not equal to number of transcripts: " << tnames.size();
            jointLog->error(errfmt.str());
            jointLog->flush();
            exit(1);
          }

          if ( eqclasses.size() != auxs_vals.size() or
               eqclasses.size() != eqclass_counts.size() ) {
            jointLog->error("Different size of the eqclasses object");
            jointLog->flush();
            exit(1);
          }
        }

        AlignmentLibraryT<UnpairedRead> alnLib(alignmentFiles,
                                               libFmt, sopt, hasEqclasses,
                                               tnames, tefflens);

        jointLog->info("Created AlignmentLibrary object");
        jointLog->flush();

        // EQCLASS
        alnLib.equivalenceClassBuilder().populateTargets(eqclasses, auxs_vals,
                                                         eqclass_counts,
                                                         alnLib.transcripts());
        success = processEqClasses(alnLib, sopt, outputDirectory);
      } else {
        AlignmentLibraryT<UnpairedRead> alnLib(alignmentFiles, transcriptFile,
                                               libFmt, sopt);

        if (autoDetectFmt) {
          alnLib.enableAutodetect();
        }
        success = processSample<UnpairedRead>(alnLib, requiredObservations, sopt,
                                              outputDirectory);
      }
    } break;
    case ReadType::PAIRED_END: {
      if (hasEqclasses) {
        jointLog->error(" Cannot quantify eqclasses in mode \n"
                        " Please report this on github");
        std::exit(1);
      }

      AlignmentLibraryT<ReadPair> alnLib(alignmentFiles, transcriptFile, libFmt,
                                        sopt);
      if (autoDetectFmt) {
        alnLib.enableAutodetect();
      }
      success = processSample<ReadPair>(alnLib, requiredObservations, sopt,
                                        outputDirectory);
    } break;
    default:
      std::stringstream errfmt;
      errfmt << "Cannot quantify library of unknown format " << libFmt;
      jointLog->error(errfmt.str());
      jointLog->flush();
      std::exit(1);
    }

    // Make sure the quantification was successful.
    if (!success) {
      jointLog->error(
          "Quantification was un-successful.  Please check the log "
          "for information about why quantification failed. If this "
          "problem persists, please report this issue on GitHub.");
      return 1;
    }

    bfs::path estFilePath = outputDirectory / "quant.sf";

    /** If the user requested gene-level abundances, then compute those now **/
    if (vm.count("geneMap")) {
      try {
        salmon::utils::generateGeneLevelEstimates(sopt.geneMapPath,
                                                  outputDirectory);
      } catch (std::exception& e) {
        fmt::print(stderr,
                   "Error: [{}] when trying to compute gene-level "
                   "estimates. The gene-level file(s) may not exist",
                   e.what());
      }
    }

  } catch (po::error& e) {
    std::cerr << "exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "============\n";
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << "============\n";
    std::cerr << argv[0] << " alignment-quant was invoked improperly.\n";
    std::cerr << "For usage information, "
              << "try " << argv[0] << " quant --help-alignments\nExiting.\n";
    std::exit(1);
  }
  return 0;
}
