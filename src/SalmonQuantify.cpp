/**
>HEADER
    Copyright (c) 2013, 2014, 2015, 2016 Rob Patro rob.patro@cs.stonybrook.edu

    This file is part of Salmon.

    Salmon is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Salmon is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Salmon.  If not, see <http://www.gnu.org/licenses/>.
<HEADER
**/

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <iterator>
#include <map>
#include <mutex>
#include <queue>
#include <random>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// C++ string formatting library
#include "spdlog/fmt/fmt.h"

// C Includes for BWA
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

// Boost Includes
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/program_options.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/thread/thread.hpp>

// Future C++ convenience classes
#include "core/range.hpp"

// TBB Includes
#include "tbb/blocked_range.h"
#include "tbb/concurrent_queue.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/partitioner.h"
#include "tbb/task_scheduler_init.h"

// logger includes
#include "spdlog/spdlog.h"

// Cereal includes
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

#include "concurrentqueue.h"

// salmon includes
#include "ClusterForest.hpp"
#include "FastxParser.hpp"
#include "IOUtils.hpp"
#include "LibraryFormat.hpp"
#include "ReadLibrary.hpp"
#include "SalmonConfig.hpp"
#include "SalmonDefaults.hpp"
#include "SalmonExceptions.hpp"
#include "SalmonIndex.hpp"
#include "SalmonMath.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"

#include "AlignmentGroup.hpp"
#include "BiasParams.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "ForgettingMassCalculator.hpp"
#include "FragmentLengthDistribution.hpp"
#include "GZipWriter.hpp"
#include "HitManager.hpp"

#include "EffectiveLengthStats.hpp"
#include "PairAlignmentFormatter.hpp"
#include "ProgramOptionsGenerator.hpp"
#include "RapMapUtils.hpp"
#include "ReadExperiment.hpp"
#include "SACollector.hpp"
#include "SASearcher.hpp"
#include "SalmonOpts.hpp"
#include "SingleAlignmentFormatter.hpp"
#include "ksw2pp/KSW2Aligner.hpp"
#include "metro/metrohash64.h"
#include "tsl/hopscotch_map.h"
#include "SelectiveAlignmentUtils.hpp"
#include "edlib.h"

/****** QUASI MAPPING DECLARATIONS *********/
using MateStatus = rapmap::utils::MateStatus;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
/****** QUASI MAPPING DECLARATIONS  *******/

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;

constexpr uint32_t miniBatchSize{5000};

template <typename AlnT> using AlnGroupVec = std::vector<AlignmentGroup<AlnT>>;

template <typename AlnT>
using AlnGroupVecRange = core::range<typename AlnGroupVec<AlnT>::iterator>;

#define __MOODYCAMEL__
#if defined(__MOODYCAMEL__)
template <typename AlnT>
using AlnGroupQueue = moodycamel::ConcurrentQueue<AlignmentGroup<AlnT>*>;
#else
template <typename AlnT>
using AlnGroupQueue = tbb::concurrent_queue<AlignmentGroup<AlnT>*>;
#endif

//#include "LightweightAlignmentDefs.hpp"

using ReadExperimentT = ReadExperiment<EquivalenceClassBuilder<TGValue>>;

template <typename AlnT>
void processMiniBatch(ReadExperimentT& readExp, ForgettingMassCalculator& fmCalc,
                      uint64_t firstTimestepOfRound, ReadLibrary& readLib,
                      const SalmonOpts& salmonOpts,
                      AlnGroupVecRange<AlnT> batchHits,
                      std::vector<Transcript>& transcripts,
                      ClusterForest& clusterForest,
                      FragmentLengthDistribution& fragLengthDist,
                      BiasParams& observedBiasParams,
                      /**
                       * NOTE : test new el model in future
                       * EffectiveLengthStats& obsEffLens,
                       */
                      std::atomic<uint64_t>& numAssignedFragments,
                      std::default_random_engine& randEng, bool initialRound,
                      std::atomic<bool>& burnedIn, double& maxZeroFrac) {

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;
  using salmon::math::LOG_EPSILON;
  using salmon::math::LOG_ONEHALF;
  using salmon::math::logAdd;
  using salmon::math::logSub;

  const uint64_t numBurninFrags = salmonOpts.numBurninFrags;

  auto& log = salmonOpts.jointLog;
  // auto log = spdlog::get("jointLog");
  size_t numTranscripts{transcripts.size()};
  size_t localNumAssignedFragments{0};
  size_t priorNumAssignedFragments{numAssignedFragments};
  std::uniform_real_distribution<> uni(
      0.0, 1.0 + std::numeric_limits<double>::min());
  std::vector<uint64_t> libTypeCounts(LibraryFormat::maxLibTypeID() + 1);
  std::vector<uint64_t> libTypeCountsPerFrag(LibraryFormat::maxLibTypeID() + 1);
  bool hasCompatibleMapping{false};
  uint64_t numCompatibleFragments{0};

  std::vector<FragmentStartPositionDistribution>& fragStartDists =
      readExp.fragmentStartPositionDistributions();
  auto& biasModel = readExp.sequenceBiasModel();
  auto& observedGCMass = observedBiasParams.observedGCMass;
  auto& obsFwd = observedBiasParams.massFwd;
  auto& obsRC = observedBiasParams.massRC;
  auto& observedPosBiasFwd = observedBiasParams.posBiasFW;
  auto& observedPosBiasRC = observedBiasParams.posBiasRC;

  bool posBiasCorrect = salmonOpts.posBiasCorrect;
  bool gcBiasCorrect = salmonOpts.gcBiasCorrect;
  bool updateCounts = initialRound;
  double incompatPrior = salmonOpts.incompatPrior;
  bool useReadCompat = incompatPrior != salmon::math::LOG_1;
  bool useFragLengthDist{!salmonOpts.noFragLengthDist};
  bool noFragLenFactor{salmonOpts.noFragLenFactor};
  bool useRankEqClasses{salmonOpts.rankEqClasses};
  uint32_t rangeFactorization{salmonOpts.rangeFactorizationBins};
  bool noLengthCorrection{salmonOpts.noLengthCorrection};
  bool useAuxParams = ((localNumAssignedFragments + numAssignedFragments) >=
                       salmonOpts.numPreBurninFrags);

  // If we're auto detecting the library type
  auto* detector = readLib.getDetector();
  bool autoDetect = (detector != nullptr) ? detector->isActive() : false;
  // If we haven't detected yet, nothing is incompatible
  if (autoDetect) {
    incompatPrior = salmon::math::LOG_1;
  }

  uint64_t zeroProbFrags{0};

  // EQClass
  auto& eqBuilder = readExp.equivalenceClassBuilder();

  // Build reverse map from transcriptID => hit id
  using HitID = uint32_t;

  double logForgettingMass{0.0};
  uint64_t currentMinibatchTimestep{0};

  // logForgettingMass and currentMinibatchTimestep are OUT parameters!
  fmCalc.getLogMassAndTimestep(logForgettingMass, currentMinibatchTimestep);

  double startingCumulativeMass =
      fmCalc.cumulativeLogMassAt(firstTimestepOfRound);

  auto isUnexpectedOrphan = [](AlnT& aln, LibraryFormat expectedLibFormat) -> bool {
    return (expectedLibFormat.type == ReadType::PAIRED_END and
            aln.mateStatus != rapmap::utils::MateStatus::PAIRED_END_PAIRED);
  };

  int i{0};
  {
    // Iterate over each group of alignments (a group consists of all alignments
    // reported
    // for a single read).  Distribute the read's mass to the transcripts
    // where it potentially aligns.
    for (auto& alnGroup : batchHits) {
      // If we had no alignments for this read, then skip it
      if (alnGroup.size() == 0) {
        continue;
      }
      LibraryFormat expectedLibraryFormat = readLib.format();
      std::fill(libTypeCountsPerFrag.begin(), libTypeCountsPerFrag.end(), 0);

      // We start out with probability 0
      double sumOfAlignProbs{LOG_0};

      // Record whether or not this read is unique to a single transcript.
      bool transcriptUnique{true};

      auto firstTranscriptID = alnGroup.alignments().front().transcriptID();
      std::unordered_set<size_t> observedTranscripts;

      // New incompat. handling.
      /**
      // The equivalence class information for
      // compatible fragments
      std::vector<uint32_t> txpIDsCompat;
      std::vector<double> auxProbsCompat;
      std::vector<double> posProbsCompat;
      double auxDenomCompat = salmon::math::LOG_0;

      // The equivalence class information for
      // all fragments (if there is no compatible fragment)
      std::vector<uint32_t> txpIDsAll;
      std::vector<double> auxProbsAll;
      std::vector<double> posProbsAll;
      double auxDenomAll = salmon::math::LOG_0;

      std::vector<uint32_t>* txpIDsFinal = nullptr;
      std::vector<uint32_t>* txpIDsFinal = nullptr;
      std::vector<uint32_t>* txpIDsFinal = nullptr;
      double auxDenomFinal = salmon::math::LOG_0;
      **/

      std::vector<uint32_t> txpIDs;
      std::vector<double> auxProbs;
      double auxDenom = salmon::math::LOG_0;

      uint32_t numInGroup{0};
      uint32_t prevTxpID{0};

      hasCompatibleMapping = false;
      // For each alignment of this read
      for (auto& aln : alnGroup.alignments()) {

        useAuxParams = ((localNumAssignedFragments + numAssignedFragments) >=
                        salmonOpts.numPreBurninFrags);
        bool considerCondProb{burnedIn or useAuxParams};

        auto transcriptID = aln.transcriptID();
        auto& transcript = transcripts[transcriptID];
        transcriptUnique =
            transcriptUnique and (transcriptID == firstTranscriptID);

        double refLength =
            transcript.RefLength > 0 ? transcript.RefLength : 1.0;
        double coverage = aln.score();
        double logFragCov = (coverage > 0) ? std::log(coverage) : LOG_1;

        // The alignment probability is the product of a
        // transcript-level term (based on abundance and) an
        // alignment-level term.
        double logRefLength{salmon::math::LOG_0};

        if (noLengthCorrection) {
          logRefLength = 1.0;
        } else if (salmonOpts.noEffectiveLengthCorrection or !burnedIn) {
          logRefLength = std::log(static_cast<double>(transcript.RefLength));
        } else {
          logRefLength = transcript.getCachedLogEffectiveLength();
        }

        double transcriptLogCount = transcript.mass(initialRound);
        auto flen = aln.fragLength();
        // If we have a properly-paired read then use the "pedantic"
        // definition here.
        if (aln.mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED and
            aln.fwd != aln.mateIsFwd) {
          flen = aln.fragLengthPedantic(transcript.RefLength);
        }

        // If the transcript had a non-zero count (including pseudocount)
        if (std::abs(transcriptLogCount) != LOG_0) {

          // The probability of drawing a fragment of this length;
          double logFragProb = LOG_1;

          // If we are expecting a paired-end library, and this is an orphan,
          // then logFragProb should be small
          if (isUnexpectedOrphan(aln, expectedLibraryFormat)) {
            logFragProb = LOG_EPSILON;
          }

          if (flen > 0.0 and useFragLengthDist and considerCondProb) {
            size_t fl = flen;
            double lenProb = fragLengthDist.pmf(fl);
            if (burnedIn) {
              /* condition fragment length prob on txp length */
              double refLengthCM =
                  fragLengthDist.cmf(static_cast<size_t>(refLength));
              bool computeMass =
                  fl < refLength and !salmon::math::isLog0(refLengthCM);
              logFragProb = (computeMass) ? (lenProb - refLengthCM)
                                          : salmon::math::LOG_EPSILON;
              if (computeMass and refLengthCM < lenProb) {
                // Threading is hard!  It's possible that an update to the PMF
                // snuck in between when we asked to cache the CMF and when the
                // "burnedIn" variable was last seen as false.
                log->info("reference length = {}, CMF[refLen] = {}, fragLen = "
                          "{}, PMF[fragLen] = {}",
                          refLength, std::exp(refLengthCM), aln.fragLength(),
                          std::exp(lenProb));
              }
            } else if (useAuxParams) {
              logFragProb = lenProb;
            }
            // logFragProb = lenProb;
            // logFragProb =
            // fragLengthDist.pmf(static_cast<size_t>(aln.fragLength()));
          }

          // TESTING
          if (noFragLenFactor) {
            logFragProb = LOG_1;
          }

          if (autoDetect) {
            detector->addSample(aln.libFormat());
            if (detector->canGuess()) {
              detector->mostLikelyType(readLib.getFormat());
              expectedLibraryFormat = readLib.getFormat();
              incompatPrior = salmonOpts.incompatPrior;
              autoDetect = false;
            } else if (!detector->isActive()) {
              expectedLibraryFormat = readLib.getFormat();
              incompatPrior = salmonOpts.incompatPrior;
              autoDetect = false;
            }
          }

          // TODO: Maybe take the fragment length distribution into account
          // for single-end fragments?

          // The probability that the fragments align to the given strands in
          // the
          // given orientations.
          bool isCompat = salmon::utils::isCompatible(
              aln.libFormat(), expectedLibraryFormat,
              static_cast<int32_t>(aln.pos), aln.fwd, aln.mateStatus);
          double logAlignCompatProb = isCompat ? LOG_1 : incompatPrior;
          if (!isCompat and salmonOpts.ignoreIncompat) {
            aln.logProb = salmon::math::LOG_0;
            continue;
          }

          /*
          double logAlignCompatProb =
              (useReadCompat) ? (salmon::utils::logAlignFormatProb(
                        aln.libFormat(), expectedLibraryFormat,
                        static_cast<int32_t>(aln.pos), aln.fwd,
                                        aln.mateStatus,
          salmonOpts.incompatPrior)) : LOG_1;
          */
          /** New compat handling
          // True if the read is compatible with the
          // expected library type; false otherwise.
          bool compat = ignoreCompat;
          if (!compat) {
              if (aln.mateStatus ==
          rapmap::utils::MateStatus::PAIRED_END_PAIRED) {
                  compat = salmon::utils::compatibleHit(
                          expectedLibType, observedLibType);
              } else {
                  int32_t pos = static_cast<int32_t>(aln.pos);
                  compat = salmon::utils::compatibleHit(
                          expectedLibraryFormat, pos,
                          aln.fwd, aln.mateStatus);
              }
          }
          **/

          // Allow for a non-uniform fragment start position distribution

          double startPosProb{-logRefLength};
          if (aln.mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED and
              !noLengthCorrection) {
            startPosProb = (flen <= refLength) ? -std::log(refLength - flen + 1)
                                               : salmon::math::LOG_EPSILON;
            // NOTE : test new el model in future
            // if (flen <= refLength) { obsEffLens.addFragment(transcriptID,
            // (refLength - flen + 1), logForgettingMass); }
          }

          double fragStartLogNumerator{salmon::math::LOG_1};
          double fragStartLogDenominator{salmon::math::LOG_1};

          auto hitPos = aln.hitPos();

          // Increment the count of this type of read that we've seen
          ++libTypeCountsPerFrag[aln.libFormat().formatID()];
          //
          if (!hasCompatibleMapping and logAlignCompatProb == LOG_1) {
            hasCompatibleMapping = true;
          }

          // The total auxiliary probabilty is the product (sum in log-space) of
          // The start position probability
          // The fragment length probabilty
          // The mapping score (coverage) probability
          // The fragment compatibility probability
          // The bias probability
          double auxProb = logFragProb + logFragCov + logAlignCompatProb;

          aln.logProb = transcriptLogCount + auxProb + startPosProb;

          // If this alignment had a zero probability, then skip it
          if (std::abs(aln.logProb) == LOG_0) {
            continue;
          }

          sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);

          if (updateCounts and observedTranscripts.find(transcriptID) ==
                                   observedTranscripts.end()) {
            transcripts[transcriptID].addTotalCount(1);
            observedTranscripts.insert(transcriptID);
          }
          // EQCLASS
          if (transcriptID < prevTxpID) {
            std::cerr << "[ERROR] Transcript IDs are not in sorted order; "
                         "please report this bug on GitHub!\n";
          }
          prevTxpID = transcriptID;
          txpIDs.push_back(transcriptID);
          auxProbs.push_back(auxProb);
          auxDenom = salmon::math::logAdd(auxDenom, auxProb);
        } else {
          aln.logProb = LOG_0;
        }
      }

      // If this fragment has a zero probability,
      // go to the next one
      if (sumOfAlignProbs == LOG_0) {
        ++zeroProbFrags;
        continue;
      } else { // otherwise, count it as assigned
        ++localNumAssignedFragments;
        if (hasCompatibleMapping) {
          ++numCompatibleFragments;
        }
      }

      // EQCLASS
      double auxProbSum{0.0};
      for (auto& p : auxProbs) {
        p = std::exp(p - auxDenom);
        auxProbSum += p;
      }

      auto eqSize = txpIDs.size();
      if (eqSize > 0) {
        if (useRankEqClasses and eqSize > 1) {
          std::vector<int> inds(eqSize);
          std::iota(inds.begin(), inds.end(), 0);
          // Get the indices in order by conditional probability
          std::sort(inds.begin(), inds.end(),
                    [&auxProbs](int i, int j) -> bool {
                      return auxProbs[i] < auxProbs[j];
                    });
          {
            decltype(txpIDs) txpIDsNew(txpIDs.size());
            decltype(auxProbs) auxProbsNew(auxProbs.size());
            for (size_t r = 0; r < eqSize; ++r) {
              auto ind = inds[r];
              txpIDsNew[r] = txpIDs[ind];
              auxProbsNew[r] = auxProbs[ind];
            }
            std::swap(txpIDsNew, txpIDs);
            std::swap(auxProbsNew, auxProbs);
          }
        }

        if (rangeFactorization > 0) {
          int32_t txpsSize = txpIDs.size();
          int32_t rangeCount = std::sqrt(txpsSize) + rangeFactorization;

          for (int32_t i = 0; i < txpsSize; i++) {
            int32_t rangeNumber = auxProbs[i] * rangeCount;
            txpIDs.push_back(rangeNumber);
          }
        }

        TranscriptGroup tg(txpIDs);
        eqBuilder.addGroup(std::move(tg), auxProbs);
      }

      // normalize the hits
      for (auto& aln : alnGroup.alignments()) {
        if (std::abs(aln.logProb) == LOG_0) {
          continue;
        }
        // Normalize the log-probability of this alignment
        aln.logProb -= sumOfAlignProbs;
        // Get the transcript referenced in this alignment
        auto transcriptID = aln.transcriptID();
        auto& transcript = transcripts[transcriptID];

        // Add the new mass to this transcript
        double newMass = logForgettingMass + aln.logProb;
        transcript.addMass(newMass);

        // Paired-end
        if (aln.libFormat().type == ReadType::PAIRED_END) {
          // TODO: Is this right for *all* library types?
          if (aln.fwd) {
            obsFwd = salmon::math::logAdd(obsFwd, aln.logProb);
          } else {
            obsRC = salmon::math::logAdd(obsRC, aln.logProb);
          }
        } else if (aln.libFormat().type == ReadType::SINGLE_END) {
          int32_t p = (aln.pos < 0) ? 0 : aln.pos;
          if (static_cast<uint32_t>(p) >= transcript.RefLength) {
            p = transcript.RefLength - 1;
          }
          // Single-end or orphan
          if (aln.libFormat().strandedness == ReadStrandedness::S) {
            obsFwd = salmon::math::logAdd(obsFwd, aln.logProb);
          } else {
            obsRC = salmon::math::logAdd(obsRC, aln.logProb);
          }
        }

        if (posBiasCorrect) {
          auto lengthClassIndex = transcript.lengthClassIndex();
          switch (aln.mateStatus) {
          case rapmap::utils::MateStatus::PAIRED_END_PAIRED: {
            // TODO: Handle the non opposite strand case
            if (aln.fwd != aln.mateIsFwd) {
              int32_t posFW = aln.fwd ? aln.pos : aln.matePos;
              int32_t posRC = aln.fwd ? aln.matePos : aln.pos;
              posFW = posFW < 0 ? 0 : posFW;
              posFW = posFW >= static_cast<int32_t>(transcript.RefLength) ?
                               static_cast<int32_t>(transcript.RefLength) - 1
                                                    : posFW;
              posRC = posRC < 0 ? 0 : posRC;
              posRC = posRC >= static_cast<int32_t>(transcript.RefLength) ?
                               static_cast<int32_t>(transcript.RefLength) - 1
                                                    : posRC;
              observedPosBiasFwd[lengthClassIndex].addMass(
                  posFW, transcript.RefLength, aln.logProb);
              observedPosBiasRC[lengthClassIndex].addMass(
                  posRC, transcript.RefLength, aln.logProb);
            }
          } break;
          case rapmap::utils::MateStatus::PAIRED_END_LEFT:
          case rapmap::utils::MateStatus::PAIRED_END_RIGHT:
          case rapmap::utils::MateStatus::SINGLE_END: {
            int32_t pos = aln.pos;
            pos = pos < 0 ? 0 : pos;
            pos = pos >= static_cast<int32_t>(transcript.RefLength) ?
                         static_cast<int32_t>(transcript.RefLength) - 1 : pos;
            if (aln.fwd) {
              observedPosBiasFwd[lengthClassIndex].addMass(
                  pos, transcript.RefLength, aln.logProb);
            } else {
              observedPosBiasRC[lengthClassIndex].addMass(
                  pos, transcript.RefLength, aln.logProb);
            }
          } break;
          default:
            break;
          }
        }

        if (gcBiasCorrect) {
          if (aln.libFormat().type == ReadType::PAIRED_END) {
            int32_t start = std::min(aln.pos, aln.matePos);
            int32_t stop = start + aln.fragLen - 1;
            // WITH CONTEXT
            if (start >= 0 and stop < static_cast<int32_t>(transcript.RefLength)) {
              bool valid{false};
              auto desc = transcript.gcDesc(start, stop, valid);
              if (valid) {
                observedGCMass.inc(desc, aln.logProb);
              }
            }
          } else if (expectedLibraryFormat.type == ReadType::SINGLE_END) {
            // Both expected and observed should be single end here
            // For single-end reads, simply assume that every fragment
            // has a length equal to the conditional mean (given the
            // current transcript's length).
            auto cmeans = readExp.condMeans();
            auto cmean =
                static_cast<int32_t>((transcript.RefLength >= cmeans.size())
                                         ? cmeans.back()
                                         : cmeans[transcript.RefLength]);
            int32_t start = aln.fwd ? aln.pos : std::max(0, aln.pos - cmean);
            int32_t stop = start + cmean;
            // WITH CONTEXT
            if (start >= 0 and stop < static_cast<int32_t>(transcript.RefLength)) {
              bool valid{false};
              auto desc = transcript.gcDesc(start, stop, valid);
              if (valid) {
                observedGCMass.inc(desc, aln.logProb);
              }
            }
          }
        }
        double r = uni(randEng);
        if (!burnedIn and r < std::exp(aln.logProb)) {

          // Old fragment length calc: double fragLength = aln.fragLength();
          auto fragLength = aln.fragLengthPedantic(transcript.RefLength);
          if (fragLength > 0) {
            fragLengthDist.addVal(fragLength, logForgettingMass);
          }

        }
      } // end normalize

      // update the single target transcript
      if (transcriptUnique) {
        if (updateCounts) {
          transcripts[firstTranscriptID].addUniqueCount(1);
        }
        clusterForest.updateCluster(firstTranscriptID, 1.0, logForgettingMass,
                                    updateCounts);
      } else { // or the appropriate clusters
        clusterForest.mergeClusters<AlnT>(alnGroup.alignments().begin(),
                                          alnGroup.alignments().end());
        clusterForest.updateCluster(
            alnGroup.alignments().front().transcriptID(), 1.0,
            logForgettingMass, updateCounts);
      }

      for(size_t i=0; i < libTypeCounts.size(); ++i) {
        libTypeCounts[i] += (libTypeCountsPerFrag[i] > 0);
      }
    } // end read group
  }   // end timer

  if (zeroProbFrags > 0) {
    auto batchReads = batchHits.size();
    maxZeroFrac = std::max(
        maxZeroFrac, static_cast<double>(100.0 * zeroProbFrags) / batchReads);
  }

  numAssignedFragments += localNumAssignedFragments;
  if (numAssignedFragments >= numBurninFrags and !burnedIn) {
    // NOTE: only one thread should succeed here, and that
    // thread will set burnedIn to true.
    readExp.updateTranscriptLengthsAtomic(burnedIn);
    fragLengthDist.cacheCMF();
  }
  if (initialRound) {
    readLib.updateLibTypeCounts(libTypeCounts);
    readLib.updateCompatCounts(numCompatibleFragments);
  }
}

using AlnCacheMap = selective_alignment::utils::AlnCacheMap;

/// START QUASI
template <typename RapMapIndexT>
void processReadsQuasi(
    paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
    AlnGroupVec<QuasiAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    RapMapIndexT* qidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
    /**
     * NOTE : test new el model in future
     * EffectiveLengthStats& obsEffLengths,
     **/
    SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache,
    MappingStatistics& mstats,
    size_t threadID) {

  uint64_t count_fwd = 0, count_bwd = 0;
  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};
  salmon::utils::ShortFragStats shortFragStats;
  double maxZeroFrac{0.0};

  // Write unmapped reads
  fmt::MemoryWriter unmappedNames;
  bool writeUnmapped = salmonOpts.writeUnmappedNames;
  spdlog::logger* unmappedLogger =
      (writeUnmapped) ? salmonOpts.unmappedLog.get() : nullptr;

  // Write unmapped reads
  fmt::MemoryWriter orphanLinks;
  bool writeOrphanLinks = salmonOpts.writeOrphanLinks;
  spdlog::logger* orphanLinkLogger =
      (writeOrphanLinks) ? salmonOpts.orphanLinkLog.get() : nullptr;

  auto& readBiasFW =
      observedBiasParams
          .seqBiasModelFW; // readExp.readBias(salmon::utils::Direction::FORWARD);
  auto& readBiasRC =
      observedBiasParams
          .seqBiasModelRC; // readExp.readBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
  // k-mers for sequence bias context
  Mer leftMer;
  Mer rightMer;

  //auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
  size_t minK = rapmap::utils::my_mer::k();

  size_t locRead{0};
  //uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};
  uint64_t localNumAssignedFragments{0};
  bool consistentHits = salmonOpts.consistentHits;
  bool quiet = salmonOpts.quiet;

  bool tooManyHits{false};
  size_t maxNumHits{salmonOpts.maxReadOccs};
  size_t readLenLeft{0};
  size_t readLenRight{0};
  SACollector<RapMapIndexT> hitCollector(qidx);

  rapmap::utils::MappingConfig mc;
  mc.consistentHits = consistentHits;
  mc.doChaining = salmonOpts.validateMappings;
  mc.consensusFraction = (salmonOpts.consensusSlack == 0.0) ? 1.0 : (1.0 - salmonOpts.consensusSlack);
  mc.allowDovetail = salmonOpts.allowDovetail;
  if (mc.doChaining) { mc.considerMultiPos = true; }

  rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<typename RapMapIndexT::IndexType>> leftHCInfo;
  rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<typename RapMapIndexT::IndexType>> rightHCInfo;

  if (salmonOpts.fasterMapping) {
    hitCollector.enableNIP();
  } else {
    hitCollector.disableNIP();
  }
  hitCollector.setStrictCheck(true);
  if (salmonOpts.quasiCoverage > 0.0) {
    hitCollector.setCoverageRequirement(salmonOpts.quasiCoverage);
  }
  if (salmonOpts.validateMappings) {
    hitCollector.enableChainScoring();
    hitCollector.setMaxMMPExtension(salmonOpts.maxMMPExtension);
  }

  bool hardFilter = salmonOpts.hardFilter;

  SASearcher<RapMapIndexT> saSearcher(qidx);
  std::vector<QuasiAlignment> leftHits;
  std::vector<QuasiAlignment> rightHits;
  rapmap::utils::HitCounters hctr;
  salmon::utils::MappingType mapType{salmon::utils::MappingType::UNMAPPED};

  PairAlignmentFormatter<RapMapIndexT*> formatter(qidx);
  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);

  std::string rc1; rc1.reserve(64);
  std::string rc2; rc2.reserve(64);

  // TODO : further investigation of bandwidth and dropoff
  using ksw2pp::KSW2Aligner;
  using ksw2pp::KSW2Config;
  using ksw2pp::EnumToType;
  using ksw2pp::KSW2AlignmentType;
  KSW2Config config;
  config.dropoff = -1;
  config.gapo = salmonOpts.gapOpenPenalty;
  config.gape = salmonOpts.gapExtendPenalty;
  config.bandwidth = salmonOpts.dpBandwidth;
  config.flag = 0;
  config.flag |= KSW_EZ_SCORE_ONLY;
  int8_t a = static_cast<int8_t>(salmonOpts.matchScore);
  int8_t b = static_cast<int8_t>(salmonOpts.mismatchPenalty);
  KSW2Aligner aligner(static_cast<int8_t>(a), static_cast<int8_t>(b));
  aligner.config() = config;
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  bool mimicStrictBT2 = salmonOpts.mimicStrictBT2;
  bool mimicBT2 = salmonOpts.mimicBT2;
  bool noDovetail = !salmonOpts.allowDovetail;
  size_t numOrphansRescued{0};

  auto ap{selective_alignment::utils::AlignmentPolicy::DEFAULT};
  if (mimicBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2;
  } else if (mimicStrictBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2_STRICT;
  }

  size_t numMappingsDropped{0};
  size_t numFragsDropped{0};
  size_t numDecoyFrags{0};

  AlnCacheMap alnCacheLeft; alnCacheLeft.reserve(32);
  AlnCacheMap alnCacheRight; alnCacheRight.reserve(32);

  auto rg = parser->getReadGroup();
  while (parser->refill(rg)) {
    rangeSize = rg.size();

    if (rangeSize > structureVec.size()) {
      salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} "
                                 "--- this shouldn't happen.\n"
                                 "Please report this bug on GitHub",
                                 rangeSize, structureVec.size());
      salmonOpts.jointLog->flush();
      spdlog::drop_all();
      std::exit(1);
    }

    for (size_t i = 0; i < rangeSize; ++i) { // For all the reads in this batch
      auto& rp = rg[i];

      // -- start resetting local variables
      // reset the status of local variables that we'll use
      // for this read.
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();
      auto& jointHits = jointHitGroup.alignments();

      leftHits.clear();
      rightHits.clear();

      leftHCInfo.clear();
      rightHCInfo.clear();
      tooManyHits = false;

      mapType = salmon::utils::MappingType::UNMAPPED;
      // -- done resetting local varaibles

      readLenLeft = rp.first.seq.length();
      readLenRight = rp.second.seq.length();
      bool tooShortLeft = (readLenLeft < minK);
      bool tooShortRight = (readLenRight < minK);

      bool lh = tooShortLeft
        ? false : hitCollector(rp.first.seq, saSearcher, leftHCInfo);


      bool rh = tooShortRight
        ? false : hitCollector(rp.second.seq, saSearcher, rightHCInfo);

      rapmap::hit_manager::hitsToMappingsSimple(*qidx, mc,
                                                MateStatus::PAIRED_END_LEFT,
                                                leftHCInfo, leftHits);

      rapmap::hit_manager::hitsToMappingsSimple(*qidx, mc,
                                                MateStatus::PAIRED_END_RIGHT,
                                                rightHCInfo, rightHits);


      rapmap::utils::MergeResult mergeRes{rapmap::utils::MergeResult::HAD_NONE};
      // Consider a read as too short if both ends are too short
      if (tooShortLeft and tooShortRight) {
        ++shortFragStats.numTooShort;
        shortFragStats.shortest = std::min(shortFragStats.shortest,
                                           std::max(readLenLeft, readLenRight));
      } else {
        rapmap::utils::MergeResult mergeRes = rapmap::utils::mergeLeftRightHitsFuzzy(lh, rh, leftHits, rightHits,
                                                                                     jointHits,
                                                                                     mc,
                                                                                     readLenLeft, maxNumHits, tooManyHits, hctr);

        // IMPORTANT NOTE : Orphan recovery currently assumes a
        // library type where mates are on separate strands
        // so (IU, ISF, ISR).  If the library type is different
        // we should either raise a warning / error, or implement
        // library-type generic recovery.
        bool mergeStatusOK = (mergeRes == rapmap::utils::MergeResult::HAD_EMPTY_INTERSECTION or
                              mergeRes == rapmap::utils::MergeResult::HAD_ONLY_LEFT or
                              mergeRes == rapmap::utils::MergeResult::HAD_ONLY_RIGHT);

        if ( mergeStatusOK and salmonOpts.recoverOrphans and !tooManyHits) {
          if (leftHits.size() + rightHits.size() > 0) {
            if (mergeRes == rapmap::utils::MergeResult::HAD_ONLY_LEFT) {
              // In this case, we've "moved" the left hit's into joint, so put them back into
              // left and make sure joint is clear.
              leftHits.swap(jointHits);
              jointHits.clear();
            } else if (mergeRes == rapmap::utils::MergeResult::HAD_ONLY_RIGHT) {
              // In this case, we've "moved" the right hit's into joint, so put them back into
              // right and make sure joint is clear.
              rightHits.swap(jointHits);
              jointHits.clear();
            }
            selective_alignment::utils::recoverOrphans(rp.first.seq,
                                                       rp.second.seq,
                                                       rc1,
                                                       rc2,
                                                       transcripts,
                                                       leftHits,
                                                       rightHits,
                                                       jointHits);
            if (!jointHits.empty()) { numOrphansRescued++; }
          }
        }

        if (initialRound) {
          upperBoundHits += (jointHits.size() > 0);
        }

        // If the read mapped to > maxReadOccs places, discard it
        if (jointHits.size() > salmonOpts.maxReadOccs) {
          jointHitGroup.clearAlignments();
        }
      }

      if (writeOrphanLinks) {
        // If we are not using strict intersection, then joint hits
        // can only be zero when either:
        // 1) there are *no* hits or
        // 2) there are hits for *both* the left and right reads, but not to the
        // same txp
        if (jointHits.size() == 0) {
          if (leftHits.size() > 0 and rightHits.size() > 0) {
            for (auto& h : leftHits) {
              orphanLinks << h.transcriptID() << ',' << h.pos << "\t";
            }
            orphanLinks << ":";
            for (auto& h : rightHits) {
              orphanLinks << h.transcriptID() << ',' << h.pos << "\t";
            }
            orphanLinks << "\n";
          }
        }
      }

      // If we have mappings, then process them.
      if (!jointHits.empty()) {
        bool isPaired = jointHits.front().mateStatus ==
                        rapmap::utils::MateStatus::PAIRED_END_PAIRED;
        if (isPaired) {
          mapType = salmon::utils::MappingType::PAIRED_MAPPED;
        }
        // If we are ignoring orphans
        if (!salmonOpts.allowOrphans) {
          // If the mappings for the current read are not properly-paired (i.e.
          // are orphans)
          // then just clear the group.
          if (!isPaired) {
            jointHitGroup.clearAlignments();
          }
        }

        bool tryAlign{salmonOpts.validateMappings};
        if (tryAlign and !jointHits.empty()) {
          alnCacheLeft.clear();
          alnCacheRight.clear();
          auto* r1 = rp.first.seq.data();
          auto* r2 = rp.second.seq.data();
          auto l1 = static_cast<int32_t>(rp.first.seq.length());
          auto l2 = static_cast<int32_t>(rp.second.seq.length());
          // We compute the reverse complements below only if we
          // need them and don't have them.
          char* r1rc = nullptr;
          char* r2rc = nullptr;
          int32_t bestScore{std::numeric_limits<int32_t>::lowest()};
          int32_t bestDecoyScore{std::numeric_limits<int32_t>::lowest()};
          std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
          std::vector<bool> decoyVec(jointHits.size(), false);
          size_t idx{0};
          double optFrac{salmonOpts.minScoreFraction};
          int32_t maxLeftScore{a * static_cast<int32_t>(rp.first.seq.length())};
          int32_t maxRightScore{a * static_cast<int32_t>(rp.second.seq.length())};
          bool multiMapping{jointHits.size() > 1};

          for (auto& h : jointHits) {
            int32_t score{std::numeric_limits<int32_t>::min()};
            auto& t = transcripts[h.tid];
            bool isDecoy = t.isDecoy();
            char* tseq = const_cast<char*>(t.Sequence());
            const int32_t tlen = static_cast<int32_t>(t.RefLength);
            const uint32_t buf{20};

            if (h.mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED) {
              if (!h.fwd and !r1rc) {
                rapmap::utils::reverseRead(rp.first.seq, rc1);
                r1rc = const_cast<char*>(rc1.data());
              }
              if (!h.mateIsFwd and !r2rc) {
                rapmap::utils::reverseRead(rp.second.seq, rc2);
                r2rc = const_cast<char*>(rc2.data());
              }
              auto* r1ptr = h.fwd ? r1 : r1rc;
              auto* r2ptr = h.mateIsFwd ? r2 : r2rc;

              int32_t s1 =
                selective_alignment::utils::getAlnScore(aligner, ez, h.pos, r1ptr, l1, tseq, tlen, a, b, maxLeftScore, h.chainStatus.getLeft(),
                            multiMapping, ap, buf, alnCacheLeft);

              int32_t s2 =
                selective_alignment::utils::getAlnScore(aligner, ez, h.matePos, r2ptr, l2, tseq, tlen, a, b, maxRightScore, h.chainStatus.getRight(),
                            multiMapping, ap, buf, alnCacheRight);

              // throw away dovetailed reads
              if (h.fwd != h.mateIsFwd and noDovetail) {
                if (h.fwd and (h.pos > h.matePos)) {
                  s1 = std::numeric_limits<int32_t>::min();
                  s2 = std::numeric_limits<int32_t>::min();
                } else if (h.mateIsFwd and (h.matePos > h.pos)) {
                  s1 = std::numeric_limits<int32_t>::min();
                  s2 = std::numeric_limits<int32_t>::min();
                }
              }

              if ((s1 < (optFrac * maxLeftScore)) or (s2 < (optFrac * maxRightScore))) {
                score = std::numeric_limits<decltype(score)>::min();
              } else {
                score = s1 + s2;
              }
            } else if (h.mateStatus == rapmap::utils::MateStatus::PAIRED_END_LEFT) {
              if (!h.fwd and !r1rc) {
                rapmap::utils::reverseRead(rp.first.seq, rc1);
                r1rc = const_cast<char*>(rc1.data());
              }
              auto* rptr = h.fwd ? r1 : r1rc;

              int32_t s =
                selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l1, tseq, tlen, a, b, maxLeftScore, h.chainStatus.getLeft(),
                            multiMapping, ap, buf, alnCacheLeft);
              if (s < (optFrac * maxLeftScore)) {
                score = std::numeric_limits<decltype(score)>::min();
              } else {
                score = s;
              }
            } else if (h.mateStatus == rapmap::utils::MateStatus::PAIRED_END_RIGHT) {
              if (!h.fwd and !r2rc) {
                rapmap::utils::reverseRead(rp.second.seq, rc2);
                r2rc = const_cast<char*>(rc2.data());
              }
              auto* rptr = h.fwd ? r2 : r2rc;

              int32_t s =
                selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l2, tseq, tlen, a, b, maxRightScore, h.chainStatus.getRight(),
                            multiMapping, ap, buf, alnCacheRight);
              if (s < (optFrac * maxRightScore)) {
                score = std::numeric_limits<decltype(score)>::min();
              } else {
                score = s;
              }
            }

            bestScore = (!isDecoy and (score > bestScore)) ? score : bestScore;
            bestDecoyScore = (isDecoy and (score > bestDecoyScore)) ? score : bestDecoyScore;
            scores[idx] = score;
            decoyVec[idx] = isDecoy;
            h.score(score);
            ++idx;
          }

          uint32_t ctr{0};
          bool bestHitDecoy = (bestScore < bestDecoyScore);
          if (bestScore > std::numeric_limits<int32_t>::min() and !bestHitDecoy) {
            // Note --- with soft filtering, only those hits that are given the minimum possible
            // score are filtered out.
            jointHits.erase(
                            std::remove_if(jointHits.begin(), jointHits.end(),
                                           [&ctr, &scores, &decoyVec, &numMappingsDropped, bestScore, hardFilter] (const QuasiAlignment& qa) -> bool {
                                             // if soft filtering, we only drop things with an invalid score
                                             // if hard filtering, we drop everything with a sub-optimal score.
                                             bool rem = decoyVec[ctr] ? true : (hardFilter ? (scores[ctr] < bestScore) :
                                                                                (scores[ctr] == std::numeric_limits<int32_t>::min()));
                                             ++ctr;
                                             numMappingsDropped += rem ? 1 : 0;
                                             return rem;
                                           }),
                            jointHits.end()
                            );

            double bestScoreD = static_cast<double>(bestScore);
            std::for_each(jointHits.begin(), jointHits.end(),
                          [bestScoreD, writeQuasimappings, hardFilter](QuasiAlignment& qa) -> void {
                            if (writeQuasimappings) { qa.alnScore(static_cast<int32_t>(qa.score())); }
                            double v = bestScoreD - qa.score();
                            qa.score( (hardFilter ? -1.0 : std::exp(-v)) );
                          });
          } else {
            numDecoyFrags += bestHitDecoy ? 1 : 0;
            ++numFragsDropped;
            jointHitGroup.clearAlignments();
          }
        } else if (isPaired and noDovetail) {
          jointHits.erase(
                          std::remove_if(jointHits.begin(), jointHits.end(),
                                         [](const QuasiAlignment& h) -> bool {
                                           if (h.fwd != h.mateIsFwd) {
                                             if (h.fwd and (h.pos > h.matePos)) {
                                               return true;
                                             } else if (h.mateIsFwd and (h.matePos > h.pos)) {
                                               return true;
                                             }
                                           }
                                           return false;
                                         }),
                          jointHits.end());
        }

        bool needBiasSample = salmonOpts.biasCorrect;

        std::uniform_int_distribution<> dis(0, jointHits.size());
        // Randomly select a hit from which to draw the bias sample.
        int32_t hitSamp{dis(eng)};
        int32_t hn{0};

        for (auto& h : jointHits) {

          // ---- Collect bias samples ------ //

          // If bias correction is turned on, and we haven't sampled a mapping
          // for this read yet, and we haven't collected the required number of
          // samples overall.
          if (needBiasSample and salmonOpts.numBiasSamples > 0 and isPaired and
              hn == hitSamp) {
            auto& t = transcripts[h.tid];

            // The "start" position is the leftmost position if
            // map to the forward strand, and the leftmost
            // position + the read length if we map to the reverse complement.

            // read 1
            int32_t pos1 = static_cast<int32_t>(h.pos);
            auto dir1 = salmon::utils::boolToDirection(h.fwd);
            int32_t startPos1 = h.fwd ? pos1 : (pos1 + h.readLen - 1);

            // read 2
            int32_t pos2 = static_cast<int32_t>(h.matePos);
            auto dir2 = salmon::utils::boolToDirection(h.mateIsFwd);
            int32_t startPos2 = h.mateIsFwd ? pos2 : (pos2 + h.mateLen - 1);

            bool success = false;

            if ((dir1 != dir2) and // Shouldn't be from the same strand
                (startPos1 > 0 and startPos1 < static_cast<int32_t>(t.RefLength)) and
                (startPos2 > 0 and startPos2 < static_cast<int32_t>(t.RefLength))) {

              const char* txpStart = t.Sequence();
              const char* txpEnd = txpStart + t.RefLength;

              const char* readStart1 = txpStart + startPos1;
              auto& readBias1 = (h.fwd) ? readBiasFW : readBiasRC;

              const char* readStart2 = txpStart + startPos2;
              auto& readBias2 = (h.mateIsFwd) ? readBiasFW : readBiasRC;

              int32_t fwPre = readBias1.contextBefore(!h.fwd);
              int32_t fwPost = readBias1.contextAfter(!h.fwd);

              int32_t rcPre = readBias2.contextBefore(!h.mateIsFwd);
              int32_t rcPost = readBias2.contextAfter(!h.mateIsFwd);

              bool read1RC = !h.fwd;
              bool read2RC = !h.mateIsFwd;

              if ((startPos1 >= readBias1.contextBefore(read1RC) and
                   startPos1 + readBias1.contextAfter(read1RC) <
                   static_cast<int32_t>(t.RefLength)) and
                  (startPos2 >= readBias2.contextBefore(read2RC) and
                   startPos2 + readBias2.contextAfter(read2RC) < static_cast<int32_t>(t.RefLength))) {

                int32_t fwPos = (h.fwd) ? startPos1 : startPos2;
                int32_t rcPos = (h.fwd) ? startPos2 : startPos1;
                if (fwPos < rcPos) {
                  leftMer.from_chars(txpStart + startPos1 -
                                     readBias1.contextBefore(read1RC));
                  rightMer.from_chars(txpStart + startPos2 -
                                      readBias2.contextBefore(read2RC));
                  if (read1RC) {
                    leftMer.reverse_complement();
                  } else {
                    rightMer.reverse_complement();
                  }

                  success = readBias1.addSequence(leftMer, 1.0);
                  success = readBias2.addSequence(rightMer, 1.0);
                }
              }

              if (success) {
                salmonOpts.numBiasSamples -= 1;
                needBiasSample = false;
              }
            }
          }
          // ---- Collect bias samples ------ //
          ++hn;

          switch (h.mateStatus) {
          case MateStatus::PAIRED_END_LEFT: {
            h.format = salmon::utils::hitType(h.pos, h.fwd);
          } break;
          case MateStatus::PAIRED_END_RIGHT: {
            // we pass in !h.fwd here because the right read
            // will have the opposite orientation from its mate.
            // NOTE : We will try recording what the mapped fragment
            // actually is, not to infer what it's mate should be.
            h.format = salmon::utils::hitType(h.pos, h.fwd);
          } break;
          case MateStatus::PAIRED_END_PAIRED: {
            uint32_t end1Pos = (h.fwd) ? h.pos : h.pos + h.readLen;
            uint32_t end2Pos =
                (h.mateIsFwd) ? h.matePos : h.matePos + h.mateLen;
            bool canDovetail = false;
            h.format =
                salmon::utils::hitType(end1Pos, h.fwd, h.readLen, end2Pos,
                                       h.mateIsFwd, h.mateLen, canDovetail);
          } break;
          case MateStatus::SINGLE_END: {
            // do nothing
          } break;
          default:
            break;
          }
        }

        if (writeQuasimappings) {
          rapmap::utils::writeAlignmentsToStream(rp, formatter, hctr, jointHits,
                                                 sstream);
        }

      } else {
        // This read was completely unmapped.
        mapType = salmon::utils::MappingType::UNMAPPED;
      }

      if (writeUnmapped and
          mapType != salmon::utils::MappingType::PAIRED_MAPPED) {
        // If we have no mappings --- then there's nothing to do
        // unless we're outputting names for un-mapped reads
        unmappedNames << rp.first.name << ' ' << salmon::utils::str(mapType)
                      << '\n';
      }

      validHits += jointHits.size();
      localNumAssignedFragments += (jointHits.size() > 0);
      locRead++;
      ++numObservedFragments;
      if (!quiet and numObservedFragments % 500000 == 0) {
        iomutex.lock();
        const char RESET_COLOR[] = "\x1b[0m";
        char green[] = "\x1b[30m";
        green[3] = '0' + static_cast<char>(fmt::GREEN);
        char red[] = "\x1b[30m";
        red[3] = '0' + static_cast<char>(fmt::RED);
        if (initialRound) {
          fmt::print(stderr, "\033[A\r\r{}processed{} {:n} {}fragments{}\n",
                     green, red, numObservedFragments, green, RESET_COLOR);
          fmt::print(stderr, "hits: {:n}, hits per frag:  {}", validHits,
                     validHits / static_cast<float>(prevObservedFrags));
        } else {
          fmt::print(stderr, "\r\r{}processed{} {:n} {}fragments{}", green, red,
                     numObservedFragments, green, RESET_COLOR);
        }
        iomutex.unlock();
      }

    } // end for i < j->nb_filled

    if (writeUnmapped) {
      std::string outStr(unmappedNames.str());
      // Get rid of last newline
      if (!outStr.empty()) {
        outStr.pop_back();
        unmappedLogger->info(std::move(outStr));
      }
      unmappedNames.clear();
    }

    if (writeQuasimappings) {
      std::string outStr(sstream.str());
      // Get rid of last newline
      if (!outStr.empty()) {
        outStr.pop_back();
        qmLog->info(std::move(outStr));
      }
      sstream.clear();
    }

    if (writeOrphanLinks) {
      std::string outStr(orphanLinks.str());
      // Get rid of last newline
      if (!outStr.empty()) {
        outStr.pop_back();
        orphanLinkLogger->info(std::move(outStr));
      }
      orphanLinks.clear();
    }

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<QuasiAlignment> hitLists = {structureVec.begin(), structureVec.begin()+rangeSize};

    /*
    if (cachingMappings) {
      salmon::utils::cacheMappings(hitLists);
    }
    */

    processMiniBatch<QuasiAlignment>(
        readExp, fmCalc, firstTimestepOfRound, rl, salmonOpts, hitLists,
        transcripts, clusterForest, fragLengthDist, observedBiasParams,
        /**
         * NOTE : test new el model in future
         * obsEffLengths,
         */
        numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac);
  }

  if (maxZeroFrac > 0.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of "
                              "{0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }

  mstats.numOrphansRescued += numOrphansRescued;
  mstats.numMappingsFiltered += numMappingsDropped;
  mstats.numFragmentsFiltered += numFragsDropped;
  mstats.numDovetails += hctr.numDovetails;
  mstats.numDecoyFragments += numDecoyFrags;
  //salmonOpts.jointLog->info("Number of orphans rescued in this thread {}",
  //                          numOrphansRescued);

  //salmonOpts.jointLog->info("Score filtering dropped {} total mappings", numDropped);
  readExp.updateShortFrags(shortFragStats);
}

// SINGLE END

// To use the parser in the following, we get ReadGroups until none is
// available.
template <typename RapMapIndexT>
void processReadsQuasi(
    single_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
    AlnGroupVec<QuasiAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    RapMapIndexT* qidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
    /**
     * NOTE : test new el model in future
     * EffectiveLengthStats& obsEffLengths,
     **/
    SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache,
    MappingStatistics& mstats,
    size_t threadID) {

  uint64_t count_fwd = 0, count_bwd = 0;
  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};
  salmon::utils::ShortFragStats shortFragStats;
  bool tooShort{false};
  double maxZeroFrac{0.0};

  // Write unmapped reads
  fmt::MemoryWriter unmappedNames;
  bool writeUnmapped = salmonOpts.writeUnmappedNames;
  spdlog::logger* unmappedLogger =
      (writeUnmapped) ? salmonOpts.unmappedLog.get() : nullptr;

  auto& readBiasFW = observedBiasParams.seqBiasModelFW;
  auto& readBiasRC = observedBiasParams.seqBiasModelRC;
  Mer context;

  const char* txomeStr = qidx->seq.c_str();

  //auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
  size_t minK = rapmap::utils::my_mer::k();

  size_t locRead{0};
  //uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};

  //bool tooManyHits{false};
  size_t readLen{0};
  size_t maxNumHits{salmonOpts.maxReadOccs};
  bool consistentHits{salmonOpts.consistentHits};
  bool quiet{salmonOpts.quiet};

  rapmap::utils::MappingConfig mc;
  mc.consistentHits = consistentHits;
  mc.doChaining = salmonOpts.validateMappings;
  mc.consensusFraction = (salmonOpts.consensusSlack == 0.0) ? 1.0 : (1.0 - salmonOpts.consensusSlack);
  mc.allowDovetail = salmonOpts.allowDovetail;
  if (mc.doChaining) { mc.considerMultiPos = true; }
  
  rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<typename RapMapIndexT::IndexType>> hcInfo;

  SACollector<RapMapIndexT> hitCollector(qidx);
  if (salmonOpts.fasterMapping) {
    hitCollector.enableNIP();
  } else {
    hitCollector.disableNIP();
  }

  hitCollector.setStrictCheck(true);
  if (salmonOpts.quasiCoverage > 0.0) {
    hitCollector.setCoverageRequirement(salmonOpts.quasiCoverage);
  }

  if (salmonOpts.validateMappings) {
    hitCollector.enableChainScoring();
    hitCollector.setMaxMMPExtension(salmonOpts.maxMMPExtension);
  }

  bool hardFilter = salmonOpts.hardFilter;

  /**
   * Setup related to mapping parameters
   **/
  SASearcher<RapMapIndexT> saSearcher(qidx);
  rapmap::utils::HitCounters hctr;

  SingleAlignmentFormatter<RapMapIndexT*> formatter(qidx);
  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);

  std::string rc1; rc1.reserve(300);

  // TODO : further investigation of bandwidth and dropoff
  using ksw2pp::KSW2Aligner;
  using ksw2pp::KSW2Config;
  using ksw2pp::EnumToType;
  using ksw2pp::KSW2AlignmentType;
  KSW2Config config;
  config.dropoff = -1;
  config.gapo = salmonOpts.gapOpenPenalty;
  config.gape = salmonOpts.gapExtendPenalty;
  config.bandwidth = salmonOpts.dpBandwidth;
  config.flag = 0;
  config.flag |= KSW_EZ_SCORE_ONLY;

  int8_t a = static_cast<int8_t>(salmonOpts.matchScore);
  int8_t b = static_cast<int8_t>(salmonOpts.mismatchPenalty);
  KSW2Aligner aligner(static_cast<int8_t>(a), static_cast<int8_t>(b));
  aligner.config() = config;
  ksw_extz_t ez;
  memset(&ez, 0, sizeof(ksw_extz_t));
  bool mimicStrictBT2 = salmonOpts.mimicStrictBT2;
  bool mimicBT2 = salmonOpts.mimicBT2;
  bool noDovetail = !salmonOpts.allowDovetail;

  auto ap{selective_alignment::utils::AlignmentPolicy::DEFAULT};
  if (mimicBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2;
  } else if (mimicStrictBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2_STRICT;
  }
  size_t numMappingsDropped{0};
  size_t numFragsDropped{0};
  size_t numDecoyFrags{0};

  //std::vector<salmon::mapping::CacheEntry> alnCache; alnCache.reserve(15);
  AlnCacheMap alnCache; alnCache.reserve(16);

  auto rg = parser->getReadGroup();
  while (parser->refill(rg)) {
    rangeSize = rg.size();
    if (rangeSize > structureVec.size()) {
      salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} "
                                 "--- this shouldn't happen.\n"
                                 "Please report this bug on GitHub",
                                 rangeSize, structureVec.size());
      salmonOpts.jointLog->flush();
      spdlog::drop_all();
      std::exit(1);
    }

    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      auto& rp = rg[i];
      readLen = rp.seq.length();
      tooShort = (readLen < minK);
      //tooManyHits = false;
      //localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      auto& jointHits = jointHitGroup.alignments();
      jointHitGroup.clearAlignments();
      hcInfo.clear();

      bool lh = tooShort ? false
        : hitCollector(rp.seq, saSearcher, hcInfo);

      rapmap::hit_manager::hitsToMappingsSimple(*qidx, mc,
                                                MateStatus::SINGLE_END,
                                                hcInfo, jointHits);

      // If the fragment was too short, record it
      if (tooShort) {
        ++shortFragStats.numTooShort;
        shortFragStats.shortest = std::min(shortFragStats.shortest, readLen);
      }

      if (initialRound) {
        upperBoundHits += (jointHits.size() > 0);
      }

      // If the read mapped to > maxReadOccs places, discard it
      if (jointHits.size() > salmonOpts.maxReadOccs) {
        jointHitGroup.clearAlignments();
      }

        bool tryAlign{salmonOpts.validateMappings};
        if (tryAlign) {
          alnCache.clear();
          auto* r1 = rp.seq.data();
          auto l1 = static_cast<int32_t>(rp.seq.length());

          char* r1rc = nullptr;
          int32_t bestScore{std::numeric_limits<int32_t>::lowest()};
          int32_t bestDecoyScore{std::numeric_limits<int32_t>::lowest()};
          std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
          std::vector<bool> decoyVec(jointHits.size(), false);
          size_t idx{0};
          double optFrac{salmonOpts.minScoreFraction};
          int32_t maxReadScore = a * rp.seq.length();
          bool multiMapping{jointHits.size() > 1};

          for (auto& h : jointHits) {
            int32_t score{std::numeric_limits<int32_t>::min()};
            auto& t = transcripts[h.tid];
            bool isDecoy = t.isDecoy();
            char* tseq = const_cast<char*>(t.Sequence());
            const int32_t tlen = static_cast<int32_t>(t.RefLength);
            const uint32_t buf{20};

            // compute the reverse complement only if we
            // need it and don't have it
            if (!h.fwd and !r1rc) {
              rapmap::utils::reverseRead(rp.seq, rc1);
              // we will not break the const promise
              r1rc = const_cast<char*>(rc1.data());
            }

            auto* rptr = h.fwd ? r1 : r1rc;
            int32_t s =
              selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l1, tseq, tlen, a, b, maxReadScore, h.chainStatus.getLeft(),
                          multiMapping, ap, buf, alnCache);
            if (s < (optFrac * maxReadScore)) {
              score = std::numeric_limits<decltype(score)>::min();
            } else {
              score = s;
            }

            bestScore = (!isDecoy and (score > bestScore)) ? score : bestScore;
            bestDecoyScore = (isDecoy and (score > bestDecoyScore)) ? score : bestDecoyScore;
            scores[idx] = score;
            decoyVec[idx] = isDecoy;
            h.score(score);
            ++idx;
          }

          uint32_t ctr{0};
          bool bestHitDecoy = (bestScore < bestDecoyScore);
          if (bestScore > std::numeric_limits<int32_t>::min() and !bestHitDecoy) {
            // Note --- with soft filtering, only those hits that are given the minimum possible
            // score are filtered out.
            jointHits.erase(
                            std::remove_if(jointHits.begin(), jointHits.end(),
                                           [&ctr, &scores, &decoyVec, &numMappingsDropped, bestScore, hardFilter] (const QuasiAlignment& qa) -> bool {
                                             // if soft filtering, we only drop things with an invalid score
                                             // if hard filtering, we drop everything with a sub-optimal score.
                                             bool rem = decoyVec[ctr] ? true : (hardFilter ? (scores[ctr] < bestScore) :
                                                                                (scores[ctr] == std::numeric_limits<int32_t>::min()));
                                             ++ctr;
                                             numMappingsDropped += rem ? 1 : 0;
                                             return rem;
                                           }),
                            jointHits.end()
                            );
            // for soft filter
            double bestScoreD = static_cast<double>(bestScore);
            std::for_each(jointHits.begin(), jointHits.end(),
                          [bestScoreD, writeQuasimappings, hardFilter](QuasiAlignment& qa) -> void {
                            if (writeQuasimappings) { qa.alnScore(static_cast<int32_t>(qa.score())); }
                            double v = bestScoreD - qa.score();
                            qa.score( (hardFilter ? -1.0 : std::exp(-v)) );
                          });
          } else {
            numDecoyFrags += bestHitDecoy ? 1 : 0;
            ++numFragsDropped;
            jointHitGroup.clearAlignments();
          }
        }


      bool needBiasSample = salmonOpts.biasCorrect;

      for (auto& h : jointHits) {

        // ---- Collect bias samples ------ //
        int32_t pos = static_cast<int32_t>(h.pos);
        auto dir = salmon::utils::boolToDirection(h.fwd);

        // If bias correction is turned on, and we haven't sampled a mapping
        // for this read yet, and we haven't collected the required number of
        // samples overall.
        if (needBiasSample and salmonOpts.numBiasSamples > 0) {
          // the "start" position is the leftmost position if
          // we hit the forward strand, and the leftmost
          // position + the read length if we hit the reverse complement
          int32_t startPos = h.fwd ? pos : pos + h.readLen;

          auto& t = transcripts[h.tid];
          if (startPos > 0 and startPos < static_cast<int32_t>(t.RefLength)) {
            auto& readBias = (h.fwd) ? readBiasFW : readBiasRC;
            const char* txpStart = t.Sequence();
            const char* readStart = txpStart + startPos;
            const char* txpEnd = txpStart + t.RefLength;

            bool success{false};
            // If the context exists around the read, add it to the observed
            // read start sequences.
            if (startPos >= readBias.contextBefore(!h.fwd) and
                startPos + readBias.contextAfter(!h.fwd) < static_cast<int32_t>(t.RefLength)) {
              context.from_chars(txpStart + startPos -
                                 readBias.contextBefore(!h.fwd));
              if (!h.fwd) {
                context.reverse_complement();
              }
              success = readBias.addSequence(context, 1.0);
            }

            if (success) {
              salmonOpts.numBiasSamples -= 1;
              needBiasSample = false;
            }
          }
        }
        // ---- Collect bias samples ------ //

        switch (h.mateStatus) {
        case MateStatus::SINGLE_END: {
          h.format = salmon::utils::hitType(h.pos, h.fwd);
        } break;
        default:
          break;
        }
      }

      if (writeQuasimappings) {
        rapmap::utils::writeAlignmentsToStream(rp, formatter, hctr, jointHits,
                                               sstream);
      }

      if (writeUnmapped and jointHits.empty()) {
        // If we have no mappings --- then there's nothing to do
        // unless we're outputting names for un-mapped reads
        unmappedNames << rp.name << " u\n";
      }

      validHits += jointHits.size();
      locRead++;
      ++numObservedFragments;
      if (!quiet and numObservedFragments % 500000 == 0) {
        iomutex.lock();
        const char RESET_COLOR[] = "\x1b[0m";
        char green[] = "\x1b[30m";
        green[3] = '0' + static_cast<char>(fmt::GREEN);
        char red[] = "\x1b[30m";
        red[3] = '0' + static_cast<char>(fmt::RED);
        if (initialRound) {
          fmt::print(stderr, "\033[A\r\r{}processed{} {:n} {}fragments{}\n",
                     green, red, numObservedFragments, green, RESET_COLOR);
          fmt::print(stderr, "hits: {:n}; hits per frag:  {}", validHits,
                     validHits / static_cast<float>(prevObservedFrags));
        } else {
          fmt::print(stderr, "\r\r{}processed{} {:n} {}fragments{}", green, red,
                     numObservedFragments, green, RESET_COLOR);
        }
        iomutex.unlock();
      }

    } // end for i < j->nb_filled

    if (writeUnmapped) {
      std::string outStr(unmappedNames.str());
      // Get rid of last newline
      if (!outStr.empty()) {
        outStr.pop_back();
        unmappedLogger->info(std::move(outStr));
      }
      unmappedNames.clear();
    }

    if (writeQuasimappings) {
      std::string outStr(sstream.str());
      // Get rid of last newline
      if (!outStr.empty()) {
        outStr.pop_back();
        qmLog->info(std::move(outStr));
      }
      sstream.clear();
    }

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<QuasiAlignment> hitLists = {structureVec.begin(), structureVec.begin()+rangeSize};
      /*boost::make_iterator_range(
        structurevec.begin(), structurevec.begin() + rangesize);*/
    processMiniBatch<QuasiAlignment>(
        readExp, fmCalc, firstTimestepOfRound, rl, salmonOpts, hitLists,
        transcripts, clusterForest, fragLengthDist, observedBiasParams,
        /**
         * NOTE : test new el model in future
         * obsEffLengths,
         **/
        numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac);
  }
  readExp.updateShortFrags(shortFragStats);

  if (maxZeroFrac > 0.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of "
                              "{0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }
  mstats.numMappingsFiltered += numMappingsDropped;
  mstats.numFragmentsFiltered += numFragsDropped;
  mstats.numDecoyFragments += numDecoyFrags;
}

/// DONE QUASI


template <typename AlnT>
void processReadLibrary(
    ReadExperimentT& readExp, ReadLibrary& rl, SalmonIndex* sidx,
    std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
    std::atomic<uint64_t>&
        numObservedFragments, // total number of reads we've looked at
    std::atomic<uint64_t>&
        numAssignedFragments,              // total number of assigned reads
    std::atomic<uint64_t>& upperBoundHits, // upper bound on # of mapped frags
    bool initialRound, std::atomic<bool>& burnedIn,
    ForgettingMassCalculator& fmCalc,
    FragmentLengthDistribution& fragLengthDist, 
    SalmonOpts& salmonOpts, double coverageThresh, bool greedyChain,
    std::mutex& iomutex, size_t numThreads,
    std::vector<AlnGroupVec<AlnT>>& structureVec, volatile bool& writeToCache, MappingStatistics& mstats) {

  std::vector<std::thread> threads;

  std::atomic<uint64_t> numValidHits{0};
  rl.checkValid();

  auto indexType = sidx->indexType();

  // Catch any exceptions that might be thrown while processing the reads
  // These two deleters are highly redundant (identical in content, but have
  // different argument types). This will be resolved by generic lambdas as soon
  // as we can rely on c++14 (waiting on bioconda).
  /** C++14 version  **/
  auto parserPtrDeleter = [&salmonOpts](auto* p) -> void {
    try {
      p->stop();
    } catch (const std::exception& e) {
      salmonOpts.jointLog->error("\n\n");
      salmonOpts.jointLog->error("Processing reads : {}", e.what());
      salmonOpts.jointLog->flush();
      spdlog::drop_all();
      std::exit(-1);
    }
    delete p;
  };
  /** C++14 version **/
  std::unique_ptr<paired_parser, decltype(parserPtrDeleter)> pairedParserPtr(
                                                                             nullptr, parserPtrDeleter);
  std::unique_ptr<single_parser, decltype(parserPtrDeleter)> singleParserPtr(
                                                                             nullptr, parserPtrDeleter);

  /** sequence-specific and GC-fragment bias vectors --- each thread gets it's
   * own **/
  std::vector<BiasParams> observedBiasParams(
      numThreads, BiasParams(salmonOpts.numConditionalGCBins,
                             salmonOpts.numFragGCBins, false));

  /**
   * NOTE : test new el model in future
   * std::vector<EffectiveLengthStats> observedEffectiveLengths(numThreads,
   *EffectiveLengthStats(numTxp));
   **/
  // NOTE : When we can support C++14, we can replace the entire ProcessFunctor class above with this
  // generic lambda.
  auto processFunctor = [&](size_t i, auto* parserPtr, auto* index) {
    if (salmonOpts.qmFileName != "" and i == 0) {
      rapmap::utils::writeSAMHeader(*index, salmonOpts.qmLog);
    }
    auto threadFun = [&, i, parserPtr, index]() -> void {
      processReadsQuasi(parserPtr, readExp, rl, structureVec[i],
                        numObservedFragments, numAssignedFragments, numValidHits,
                        upperBoundHits, index, transcripts,
                        fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                        salmonOpts, coverageThresh, iomutex, initialRound,
                        burnedIn, writeToCache, mstats, i);
    };
    threads.emplace_back(threadFun);
  };

  // If the read library is paired-end
  // ------ Paired-end --------
  bool isPairedEnd = rl.format().type == ReadType::PAIRED_END;
  bool isSingleEnd = rl.format().type == ReadType::SINGLE_END;

  if (isPairedEnd) {

    if (rl.mates1().size() != rl.mates2().size()) {
      salmonOpts.jointLog->error("The number of provided files for "
                                 "-1 and -2 must be the same!");
      std::exit(1);
    }

    size_t numFiles = rl.mates1().size() + rl.mates2().size();
    uint32_t numParsingThreads{1};
    // HACK!
    if (rl.mates1().size() > 1 and numThreads > 8) {
      numParsingThreads = 2;
    }
    pairedParserPtr.reset(new paired_parser(rl.mates1(), rl.mates2(),
                                            numThreads, numParsingThreads,
                                            miniBatchSize));
    pairedParserPtr->start();
  } else if (isSingleEnd) {
    uint32_t numParsingThreads{1};
    // HACK!
    if (rl.unmated().size() > 1 and numThreads > 8) {
      numParsingThreads = 2;
    }
    singleParserPtr.reset(new single_parser(rl.unmated(), numThreads,
                                            numParsingThreads, miniBatchSize));
    singleParserPtr->start();
  }

  switch (indexType) {
    case SalmonIndexType::FMD: {
      fmt::MemoryWriter infostr;
      infostr << "This version of salmon does not support FMD indexing.";
      throw std::invalid_argument(infostr.str());
    } break;
    case SalmonIndexType::QUASI: {
      // True if we have a 64-bit SA index, false otherwise
      bool largeIndex = sidx->is64BitQuasi();
      bool perfectHashIndex = sidx->isPerfectHashQuasi();
      for (size_t i = 0; i < numThreads; ++i) {
        // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
        // change value before the lambda below is evaluated --- crazy!
        if (largeIndex) {
          if (perfectHashIndex) { // Perfect Hash
            if (isPairedEnd) {processFunctor(i, pairedParserPtr.get(), sidx->quasiIndexPerfectHash64());}
            else if (isSingleEnd) {processFunctor(i, singleParserPtr.get(), sidx->quasiIndexPerfectHash64());}
          } else { // Dense Hash
            if (isPairedEnd) {processFunctor(i, pairedParserPtr.get(), sidx->quasiIndex64());}
            else if (isSingleEnd) {processFunctor(i, singleParserPtr.get(), sidx->quasiIndex64());}
          }
        } else {
          if (perfectHashIndex) { // Perfect Hash
            if (isPairedEnd) { processFunctor(i, pairedParserPtr.get(), sidx->quasiIndexPerfectHash32()); }
            else if (isSingleEnd) { processFunctor(i, singleParserPtr.get(), sidx->quasiIndexPerfectHash32()); }
          } else { // Dense Hash
            if (isPairedEnd) { processFunctor(i, pairedParserPtr.get(), sidx->quasiIndex32()); }
            else if (isSingleEnd) { processFunctor(i, singleParserPtr.get(), sidx->quasiIndex32()); }
          }
        } // End spawn current thread

      } // End spawn all threads
    }   // End Quasi index
    break;
  } // end switch

  for (auto& t : threads) {
    t.join();
  }

  // At this point, if we were using decoy transcripts, we don't need them anymore and can get
  // rid of them.
  readExp.dropDecoyTranscripts();

  // If we don't have a sufficient number of assigned fragments, then
  // complain here!
  if (numAssignedFragments < salmonOpts.minRequiredFrags) {
    readExp.setNumObservedFragments(numObservedFragments);
    readExp.numAssignedFragmentsAtomic().store(numAssignedFragments);
    double mappingRate = numAssignedFragments.load() /
                          static_cast<double>(numObservedFragments.load());
    readExp.setEffectiveMappingRate(mappingRate);
    throw InsufficientAssignedFragments(numAssignedFragments.load(),
                                        salmonOpts.minRequiredFrags);
  }

  /**
    * NOTE : test new el model in future
  EffectiveLengthStats eel(numTxp);
  for (auto& els : observedEffectiveLengths) {
    eel.merge(els);
  }
  auto& transcripts = readExp.transcripts();
  for (size_t tid = 0; tid < numTxp; ++tid) {
    auto el = eel.getExpectedEffectiveLength(tid);
    auto countObs = eel.getObservedCount(tid);
    if (countObs > salmonOpts.eelCountCutoff and el >= 1.0) {
      transcripts[tid].setCachedLogEffectiveLength(std::log(el));
    }
  }
  **/

  /** GC-fragment bias **/
  // Set the global distribution based on the sum of local
  // distributions.
  double gcFracFwd{0.0};
  double globalMass{salmon::math::LOG_0};
  double globalFwdMass{salmon::math::LOG_0};
  auto& globalGCMass = readExp.observedGC();
  for (auto& gcp : observedBiasParams) {
    auto& gcm = gcp.observedGCMass;
    globalGCMass.combineCounts(gcm);

    auto& fw =
        readExp.readBiasModelObserved(salmon::utils::Direction::FORWARD);
    auto& rc = readExp.readBiasModelObserved(
        salmon::utils::Direction::REVERSE_COMPLEMENT);

    auto& fwloc = gcp.seqBiasModelFW;
    auto& rcloc = gcp.seqBiasModelRC;
    fw.combineCounts(fwloc);
    rc.combineCounts(rcloc);

    /**
      * positional biases
      **/
    auto& posBiasesFW = readExp.posBias(salmon::utils::Direction::FORWARD);
    auto& posBiasesRC =
        readExp.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
    for (size_t i = 0; i < posBiasesFW.size(); ++i) {
      posBiasesFW[i].combine(gcp.posBiasFW[i]);
      posBiasesRC[i].combine(gcp.posBiasRC[i]);
    }
    /*
            for (size_t i = 0; i < fwloc.counts.size(); ++i) {
                fw.counts[i] += fwloc.counts[i];
                rc.counts[i] += rcloc.counts[i];
            }
    */

    globalMass = salmon::math::logAdd(globalMass, gcp.massFwd);
    globalMass = salmon::math::logAdd(globalMass, gcp.massRC);
    globalFwdMass = salmon::math::logAdd(globalFwdMass, gcp.massFwd);
  }
  globalGCMass.normalize();

  if (globalMass != salmon::math::LOG_0) {
    if (globalFwdMass != salmon::math::LOG_0) {
      gcFracFwd = std::exp(globalFwdMass - globalMass);
    }
    readExp.setGCFracForward(gcFracFwd);
  }

  // finalize the positional biases
  if (salmonOpts.posBiasCorrect) {
    auto& posBiasesFW = readExp.posBias(salmon::utils::Direction::FORWARD);
    auto& posBiasesRC =
        readExp.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
    for (size_t i = 0; i < posBiasesFW.size(); ++i) {
      posBiasesFW[i].finalize();
      posBiasesRC[i].finalize();
    }
  }
  /** END GC-fragment bias **/
}

/**
 *  Quantify the targets given in the file `transcriptFile` using the
 *  reads in the given set of `readLibraries`, and write the results
 *  to the file `outputFile`.  The reads are assumed to be in the format
 *  specified by `libFmt`.
 *
 */
template <typename AlnT>
void quantifyLibrary(ReadExperimentT& experiment, bool greedyChain,
                     SalmonOpts& salmonOpts,
                     MappingStatistics& mstats,
                     double coverageThresh, uint32_t numQuantThreads) {

  bool burnedIn = (salmonOpts.numBurninFrags == 0);
  uint64_t numRequiredFragments = salmonOpts.numRequiredFragments;
  std::atomic<uint64_t> upperBoundHits{0};
  auto& refs = experiment.transcripts();
  size_t numTranscripts = refs.size();
  // The *total* number of fragments observed so far (over all passes through
  // the data).
  std::atomic<uint64_t> numObservedFragments{0};
  uint64_t prevNumObservedFragments{0};
  // The *total* number of fragments assigned so far (over all passes through
  // the data).
  std::atomic<uint64_t> totalAssignedFragments{0};
  uint64_t prevNumAssignedFragments{0};

  auto jointLog = salmonOpts.jointLog;

  ForgettingMassCalculator fmCalc(salmonOpts.forgettingFactor);
  size_t prefillSize = 1000000000 / miniBatchSize;
  fmCalc.prefill(prefillSize);

  bool initialRound{true};
  uint32_t roundNum{0};

  std::mutex ffMutex;
  std::mutex ioMutex;

  size_t numPrevObservedFragments = 0;

  size_t maxReadGroup{miniBatchSize};
  uint32_t structCacheSize = numQuantThreads * maxReadGroup * 10;

  // EQCLASS
  bool terminate{false};

  while (numObservedFragments < numRequiredFragments and !terminate) {
    prevNumObservedFragments = numObservedFragments;
    if (!initialRound) {
      bool didReset = (salmonOpts.disableMappingCache)
                          ? (experiment.reset())
                          : (experiment.softReset());

      if (!didReset) {
        std::string errmsg = fmt::format(
            "\n\n======== WARNING ========\n"
            "One of the provided read files: [{}] "
            "is not a regular file and therefore can't be read from "
            "more than once.\n\n"
            "We observed only {} mapping fragments when we wanted at least "
            "{}.\n\n"
            "Please consider re-running Salmon with these reads "
            "as a regular file!\n"
            "NOTE: If you received this warning from salmon but did not "
            "disable the mapping cache (--disableMappingCache), then there \n"
            "was some other problem. Please make sure, e.g., that you have not "
            "run out of disk space.\n"
            "==========================\n\n",
            experiment.readFilesAsString(), numObservedFragments,
            numRequiredFragments);
        jointLog->warn(errmsg);
        break;
      }

      numPrevObservedFragments = numObservedFragments;
    }

    // This structure is a vector of vectors of alignment
    // groups.  Each thread will get its own vector, so we
    // allocate these up front to save time and allow
    // reuse.
    std::vector<AlnGroupVec<AlnT>> groupVec;
    for (size_t i = 0; i < numQuantThreads; ++i) {
      groupVec.emplace_back(maxReadGroup);
    }

    bool writeToCache = !salmonOpts.disableMappingCache;
    auto processReadLibraryCallback =
        [&](ReadLibrary& rl, SalmonIndex* sidx,
            std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
            FragmentLengthDistribution& fragLengthDist,
            std::atomic<uint64_t>& numAssignedFragments, size_t numQuantThreads,
            std::atomic<bool>& burnedIn) -> void {

      processReadLibrary<AlnT>(experiment, rl, sidx, transcripts, clusterForest,
                               numObservedFragments, totalAssignedFragments,
                               upperBoundHits, initialRound, burnedIn, fmCalc,
                               fragLengthDist, salmonOpts,
                               coverageThresh, greedyChain, ioMutex,
                               numQuantThreads, groupVec, writeToCache, mstats);

      numAssignedFragments = totalAssignedFragments - prevNumAssignedFragments;
      prevNumAssignedFragments = totalAssignedFragments;
    };

    // Process all of the reads
    if (!salmonOpts.quiet) {
      fmt::print(stderr, "\n\n\n\n");
    }
    experiment.processReads(numQuantThreads, salmonOpts,
                            processReadLibraryCallback);
    experiment.setNumObservedFragments(numObservedFragments);

    // EQCLASS
    bool done = experiment.equivalenceClassBuilder().finish();
    // skip the extra online rounds
    terminate = true;

    initialRound = false;
    ++roundNum;

    if (!salmonOpts.quiet) {
      fmt::print(stderr, "\n\n\n\n");
    }
    /*
    fmt::print(stderr, "\n# observed = {} / # required = {}\n",
               numObservedFragments, numRequiredFragments);
    fmt::print(stderr, "hard # assigned = {} / # observed (this round) = {} : "
                       "upper bound assigned = {} \033[A\033[A",
               experiment.numAssignedFragments(),
               numObservedFragments - numPrevObservedFragments,
               upperBoundHits);
    */
    salmonOpts.fileLog->info(
        "\nAt end of round {}\n"
        "==================\n"
        "Observed {} total fragments ({} in most recent round)\n",
        roundNum - 1, numObservedFragments,
        numObservedFragments - numPrevObservedFragments);
  }
  if (!salmonOpts.quiet) {
    fmt::print(stderr, "\n\n\n\n");
  }

  // Report statistics about short fragments
  salmon::utils::ShortFragStats shortFragStats = experiment.getShortFragStats();
  if (shortFragStats.numTooShort > 0) {
    double tooShortFrac =
        (numObservedFragments > 0)
            ? (static_cast<double>(shortFragStats.numTooShort) /
               numObservedFragments)
            : 0.0;
    if (tooShortFrac > 0.0) {
      size_t minK = rapmap::utils::my_mer::k();
      fmt::print(stderr, "\n\n");
      salmonOpts.jointLog->warn("{}% of fragments were shorter than the k used "
                                "to build the index ({}).\n"
                                "If this fraction is too large, consider "
                                "re-building the index with a smaller k.\n"
                                "The minimum read size found was {}.\n\n",
                                tooShortFrac * 100.0, minK,
                                shortFragStats.shortest);

      // If *all* fragments were too short, then halt now
      if (shortFragStats.numTooShort == numObservedFragments) {
        salmonOpts.jointLog->error(
            "All fragments were too short to quasi-map.  I won't proceed.");
        std::exit(1);
      }
    } // end tooShortFrac > 0.0
  }

  if (salmonOpts.recoverOrphans) {
    salmonOpts.jointLog->info("Number of orphans recovered using orphan rescue : {:n}", mstats.numOrphansRescued.load());
  }
  if (salmonOpts.validateMappings) {
    salmonOpts.jointLog->info("Number of mappings discarded because of alignment score : {:n}", mstats.numMappingsFiltered.load());
    salmonOpts.jointLog->info("Number of fragments entirely discarded because of alignment score : {:n}", mstats.numFragmentsFiltered.load());
    salmonOpts.jointLog->info("Number of fragments discarded because they are best-mapped to decoys : {:n}", mstats.numDecoyFragments.load());
  }
  if (!salmonOpts.allowDovetail) {
    salmonOpts.jointLog->info("Number of fragments discarded because they have only dovetail (discordant) mappings : {:n}", mstats.numDovetails.load());
  }

  // If we didn't achieve burnin, then at least compute effective
  // lengths and mention this to the user.
  if (totalAssignedFragments < salmonOpts.numBurninFrags) {
    std::atomic<bool> dummyBool{false};
    experiment.updateTranscriptLengthsAtomic(dummyBool);

    jointLog->warn("Only {} fragments were mapped, but the number of burn-in "
                   "fragments was set to {}.\n"
                   "The effective lengths have been computed using the "
                   "observed mappings.\n",
                   totalAssignedFragments, salmonOpts.numBurninFrags);
  }

  if (numObservedFragments <= prevNumObservedFragments) {
    jointLog->warn(
        "Something seems to be wrong with the calculation "
        "of the mapping rate.  The recorded ratio is likely wrong.  Please "
        "file this as a bug report.\n");
  } else {
    double upperBoundMappingRate =
        upperBoundHits.load() /
        static_cast<double>(numObservedFragments.load());
    experiment.setNumObservedFragments(numObservedFragments -
                                       prevNumObservedFragments);
    experiment.setUpperBoundHits(upperBoundHits.load());
    double mappingRate = totalAssignedFragments.load() /
      static_cast<double>(numObservedFragments.load());
    experiment.setEffectiveMappingRate(mappingRate);
  }

  jointLog->info("Mapping rate = {}\%\n",
                 experiment.effectiveMappingRate() * 100.0);
  jointLog->info("finished quantifyLibrary()");
}

int salmonQuantify(int argc, const char* argv[]) {
  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  int32_t numBiasSamples{0};

  SalmonOpts sopt;

  sopt.numThreads = std::thread::hardware_concurrency();

  salmon::ProgramOptionsGenerator pogen;

  auto inputOpt = pogen.getMappingInputOptions(sopt);
  auto basicOpt = pogen.getBasicOptions(sopt);
  auto mapSpecOpt = pogen.getMappingSpecificOptions(sopt);
  auto advancedOpt = pogen.getAdvancedOptions(numBiasSamples, sopt);
  auto hiddenOpt = pogen.getHiddenOptions(sopt);
  auto testingOpt = pogen.getTestingOptions(sopt);
  auto deprecatedOpt = pogen.getDeprecatedOptions(sopt);

  po::options_description all("salmon quant options");
  all.add(inputOpt).add(basicOpt).add(mapSpecOpt).add(advancedOpt).add(testingOpt).add(hiddenOpt).add(deprecatedOpt);

  po::options_description visible("salmon quant options");
  visible.add(inputOpt).add(basicOpt).add(mapSpecOpt).add(advancedOpt);

  po::variables_map vm;
  try {
    auto orderedOptions =
        po::command_line_parser(argc, argv).options(all).run();

    po::store(orderedOptions, vm);

    if (vm.count("help")) {
      auto hstring = R"(
Quant
==========
Perform dual-phase, mapping-based estimation of
transcript abundance from RNA-seq reads
)";
      std::cout << hstring << std::endl;
      std::cout << visible << std::endl;
      std::exit(0);
    }

    po::notify(vm);

    // If we're supposed to be quiet, set the global logger level to >= warn
    if (sopt.quiet) {
      spdlog::set_level(spdlog::level::warn); // Set global log level to info
    }

    std::stringstream commentStream;
    commentStream << "### salmon (mapping-based) v" << salmon::version << "\n";
    commentStream << "### [ program ] => salmon \n";
    commentStream << "### [ command ] => quant \n";
    for (auto& opt : orderedOptions.options) {
      commentStream << "### [ " << opt.string_key << " ] => {";
      for (auto& val : opt.value) {
        commentStream << " " << val;
      }
      commentStream << " }\n";
    }
    std::string commentString = commentStream.str();
    if (!sopt.quiet) {
      fmt::print(stderr, "{}", commentString);
    }

    sopt.quantMode = SalmonQuantMode::MAP;
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
    bool greedyChain = true;

    jointLog->info("parsing read library format");

    // ==== Library format processing ===
    vector<ReadLibrary> readLibraries =
      salmon::utils::extractReadLibraries(orderedOptions);

    if (readLibraries.size() == 0) {
      jointLog->error(
          "Failed to successfully parse any complete read libraries. "
          "Please make sure you provided arguments properly to -1, -2 (for "
          "paired-end libraries) "
          "or -r (for single-end libraries), and that the library format "
          "option (-l) *comes before* the read libraries.");
      std::exit(1);
    }
    // ==== END: Library format processing ===

    SalmonIndexVersionInfo versionInfo;
    boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
    versionInfo.load(versionPath);
    auto idxType = versionInfo.indexType();

    MappingStatistics mstats;
    ReadExperimentT experiment(readLibraries, indexDirectory, sopt);

    // This will be the class in charge of maintaining our
    // rich equivalence classes
    experiment.equivalenceClassBuilder().setMaxResizeThreads(sopt.maxHashResizeThreads);
    experiment.equivalenceClassBuilder().start();

    auto indexType = experiment.getIndex()->indexType();

    try {
      switch (indexType) {
      case SalmonIndexType::FMD: {
        fmt::MemoryWriter infostr;
        infostr << "This version of salmon does not support FMD indexing.";
        throw std::invalid_argument(infostr.str());
      } break;
      case SalmonIndexType::QUASI: {
        // We can only do fragment GC bias correction, for the time being, with
        // paired-end reads
        if (sopt.gcBiasCorrect) {
          for (auto& rl : readLibraries) {
            if (rl.format().type != ReadType::PAIRED_END) {
              jointLog->warn(
                  "Fragment GC bias correction is currently *experimental* "
                  "in single-end libraries.  Please use this option "
                  "with caution.");
              // sopt.gcBiasCorrect = false;
            }
          }
        }

        sopt.allowOrphans = !sopt.discardOrphansQuasi;
        sopt.useQuasi = true;
        quantifyLibrary<QuasiAlignment>(experiment, greedyChain,
                                        sopt, mstats, sopt.coverageThresh, sopt.numThreads);
      } break;
      }
    } catch (const InsufficientAssignedFragments& iaf) {
      sopt.jointLog->warn(iaf.what());
      salmon::utils::writeCmdInfo(sopt, orderedOptions);
      GZipWriter gzw(outputDirectory, jointLog);
      gzw.writeEmptyAbundances(sopt, experiment);
      // Write meta-information about the run
      std::vector<std::string> errors{"insufficient_assigned_fragments"};
      sopt.runStopTime = salmon::utils::getCurrentTimeAsString();
      gzw.writeEmptyMeta(sopt, experiment, errors);
      std::exit(1);
    }

    // Write out information about the command / run
    salmon::utils::writeCmdInfo(sopt, orderedOptions);

    GZipWriter gzw(outputDirectory, jointLog);

    if (!sopt.skipQuant) {
      // Now that the streaming pass is complete, we have
      // our initial estimates, and our rich equivalence
      // classes.  Perform further optimization until
      // convergence.
      // NOTE: A side-effect of calling the optimizer is that
      // the `EffectiveLength` field of each transcript is
      // set to its final value.
      CollapsedEMOptimizer optimizer;
      jointLog->info("Starting optimizer");
      salmon::utils::normalizeAlphas(sopt, experiment);
      bool optSuccess = optimizer.optimize(experiment, sopt, 0.01, 10000);

      if (!optSuccess) {
        jointLog->error(
                        "The optimization algorithm failed. This is likely the result of "
                        "bad input (or a bug). If you cannot track down the cause, please "
                        "report this issue on GitHub.");
        return 1;
      }
      jointLog->info("Finished optimizer");

      jointLog->info("writing output \n");

      // Write the quantification results
      gzw.writeAbundances(sopt, experiment);

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
          // sampler.sampleMultipleChains(experiment, sopt, bsWriter,
          // sopt.numGibbsSamples);
          sampler.sample(experiment, sopt, bsWriter, sopt.numGibbsSamples);
        if (!sampleSuccess) {
          jointLog->error("Encountered error during Gibbs sampling.\n"
                          "This should not happen.\n"
                          "Please file a bug report on GitHub.\n");
          return 1;
        }
        jointLog->info("Finished Gibbs Sampler");
      } else if (sopt.numBootstraps > 0) {
        gzw.setSamplingPath(sopt);
        // The function we'll use as a callback to write samples
        std::function<bool(const std::vector<double>&)> bsWriter =
          [&gzw](const std::vector<double>& alphas) -> bool {
            return gzw.writeBootstrap(alphas);
          };

        jointLog->info("Starting Bootstrapping");
        bool bootstrapSuccess =
          optimizer.gatherBootstraps(experiment, sopt, bsWriter, 0.01, 10000);
        jointLog->info("Finished Bootstrapping");
        if (!bootstrapSuccess) {
          jointLog->error("Encountered error during bootstrapping.\n"
                          "This should not happen.\n"
                          "Please file a bug report on GitHub.\n");
          return 1;
        }
      }

      /** If the user requested gene-level abundances, then compute those now **/
      if (vm.count("geneMap")) {
        try {
          salmon::utils::generateGeneLevelEstimates(sopt.geneMapPath,
                                                    outputDirectory);
        } catch (std::invalid_argument& e) {
          fmt::print(stderr,
                     "Error: [{}] when trying to compute gene-level "
                     "estimates. The gene-level file(s) may not exist",
                     e.what());
        }
      }
    }

    // If we are dumping the equivalence classes, then
    // do it here.
    if (sopt.dumpEq) {
      jointLog->info("writing equivalence class counts.");
      gzw.writeEquivCounts(sopt, experiment);
      jointLog->info("done writing equivalence class counts.");
    }

    bfs::path libCountFilePath = outputDirectory / "lib_format_counts.json";
    experiment.summarizeLibraryTypeCounts(libCountFilePath);

    // Test writing out the fragment length distribution
    if (!sopt.noFragLengthDist) {
      bfs::path distFileName = sopt.paramsDirectory / "flenDist.txt";
      {
        std::unique_ptr<std::FILE, int (*)(std::FILE*)> distOut(
                                                                std::fopen(distFileName.c_str(), "w"), std::fclose);
        fmt::print(distOut.get(), "{}\n",
                   experiment.fragmentLengthDistribution()->toString());
      }
    }

    if (sopt.writeUnmappedNames) {
      auto l = sopt.unmappedLog.get();
      // If the logger was created, then flush it and
      // close the associated file.
      if (l) {
        l->flush();
        if (sopt.unmappedFile) {
          sopt.unmappedFile->close();
        }
      }
    }

    if (sopt.writeOrphanLinks) {
      auto l = sopt.orphanLinkLog.get();
      // If the logger was created, then flush it and
      // close the associated file.
      if (l) {
        l->flush();
        if (sopt.orphanLinkFile) {
          sopt.orphanLinkFile->close();
        }
      }
    }

    // if we wrote quasimappings, flush that buffer
    if (sopt.qmFileName != "") {
      sopt.qmLog->flush();
      // if we wrote to a buffer other than stdout, close
      // the file
      if (sopt.qmFileName != "-") {
        sopt.qmFile.close();
      }
    }

    sopt.runStopTime = salmon::utils::getCurrentTimeAsString();

    // Write meta-information about the run
    gzw.writeMeta(sopt, experiment, mstats);

  } catch (po::error& e) {
    std::cerr << "(mapping-based mode) Exception : [" << e.what() << "].\n";
    std::cerr << "Please be sure you are passing correct options, and that you are running in the intended mode.\n";
    std::cerr << "alignment-based mode is detected and enabled via the \'-a\' flag. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << argv[0] << " quant was invoked improperly.\n";
    std::cerr << "For usage information, try " << argv[0]
              << " quant --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}
