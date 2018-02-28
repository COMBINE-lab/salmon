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


// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"

// Boost Includes
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/program_options.hpp>
#include <boost/range/irange.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/thread/thread.hpp>

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

#include <cuckoohash_map.hh>

//alevin include
#include "AlevinOpts.hpp"
#include "AlevinUtils.hpp"
#include "SingleCellProtocols.hpp"
#include "CollapsedCellOptimizer.hpp"
#include "BarcodeGroup.hpp"
#include "Dedup.hpp"

// salmon includes
#include "ClusterForest.hpp"
#include "FastxParser.hpp"
#include "IOUtils.hpp"
#include "LibraryFormat.hpp"
#include "ReadLibrary.hpp"
#include "SalmonConfig.hpp"
#include "SalmonIndex.hpp"
#include "SalmonMath.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"

#include "AlignmentGroup.hpp"
#include "BWAUtils.hpp"
#include "BiasParams.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "ForgettingMassCalculator.hpp"
#include "FragmentLengthDistribution.hpp"
#include "GZipWriter.hpp"
#include "HitManager.hpp"
#include "KmerIntervalMap.hpp"

#include "RapMapUtils.hpp"
#include "ReadExperiment.hpp"
#include "SACollector.hpp"
#include "SASearcher.hpp"
#include "SalmonOpts.hpp"
#include "PairAlignmentFormatter.hpp"
#include "SingleAlignmentFormatter.hpp"
#include "RapMapUtils.hpp"
#include "AlevinKmer.hpp"

namespace alevin{
  extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "ksort.h"
#include "kvec.h"
#include "utils.h"
  }

  /****** QUASI MAPPING DECLARATIONS *********/
  using MateStatus = rapmap::utils::MateStatus;
  using QuasiAlignment = rapmap::utils::QuasiAlignment;
  /****** QUASI MAPPING DECLARATIONS  *******/

  using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  using TranscriptID = uint32_t;
  using TranscriptIDVector = std::vector<TranscriptID>;
  using KmerIDMap = std::vector<TranscriptIDVector>;
  using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

  constexpr uint32_t miniBatchSize{5000};

  using CellBarcodeT = uint32_t;
  using UMIBarcodeT = uint64_t;

  template <typename AlnT> using AlevinAlnGroup = AlignmentGroup<AlnT, CellBarcodeT, UMIBarcodeT>;
  template <typename AlnT> using AlnGroupVec = std::vector<AlevinAlnGroup<AlnT>>;

  template <typename AlnT>
  using AlnGroupVecRange =
    boost::iterator_range<typename AlnGroupVec<AlnT>::iterator>;

#define __MOODYCAMEL__
#if defined(__MOODYCAMEL__)
  template <typename AlnT>
  using AlnGroupQueue = moodycamel::ConcurrentQueue<AlevinAlnGroup<AlnT>*>;
#else
  template <typename AlnT>
  using AlnGroupQueue = tbb::concurrent_queue<AlevinAlnGroup<AlnT>*>;
#endif

#include "LightweightAlignmentDefs.hpp"
}

//have to create new namespace because of multiple definition
using namespace alevin;

/* ALEVIN DECLERATIONS*/
using bcEnd = BarcodeEnd;
namespace aut = alevin::utils;
using BlockedIndexRange = tbb::blocked_range<size_t>;
using ReadExperimentT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;


template <typename AlnT, typename ProtocolT>
void processMiniBatch(ReadExperimentT& readExp, ForgettingMassCalculator& fmCalc,
                              uint64_t firstTimestepOfRound, ReadLibrary& readLib,
                              const SalmonOpts& salmonOpts,
                              AlnGroupVecRange<AlnT> batchHits,
                              std::vector<Transcript>& transcripts,
                              ClusterForest& clusterForest,
                              FragmentLengthDistribution& fragLengthDist,
                              BiasParams& observedBiasParams,
                              std::atomic<uint64_t>& numAssignedFragments,
                              std::default_random_engine& randEng, bool initialRound,
                              std::atomic<bool>& burnedIn, double& maxZeroFrac,
                              AlevinOpts<ProtocolT>& alevinOpts,
                              std::mutex& iomutex){

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;
  using salmon::math::LOG_EPSILON;
  using salmon::math::LOG_ONEHALF;
  using salmon::math::logAdd;
  using salmon::math::logSub;

  const uint64_t numBurninFrags = salmonOpts.numBurninFrags;

  auto& log = salmonOpts.jointLog;
  //auto log = spdlog::get("jointLog");
  size_t numTranscripts{transcripts.size()};
  size_t localNumAssignedFragments{0};
  size_t priorNumAssignedFragments{numAssignedFragments};
  std::uniform_real_distribution<> uni(
                                       0.0, 1.0 + std::numeric_limits<double>::min());
  std::vector<uint64_t> libTypeCounts(LibraryFormat::maxLibTypeID() + 1);
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
  bool useFSPD{salmonOpts.useFSPD};
  bool useFragLengthDist{!salmonOpts.noFragLengthDist};
  bool noFragLenFactor{salmonOpts.noFragLenFactor};
  bool useRankEqClasses{salmonOpts.rankEqClasses};
  bool noLengthCorrection{salmonOpts.noLengthCorrection};
  // JAN 13
  bool useAuxParams = ((localNumAssignedFragments + numAssignedFragments) >= salmonOpts.numPreBurninFrags);

  // If we're auto detecting the library type
  auto* detector = readLib.getDetector();
  bool autoDetect = (detector != nullptr) ? detector->isActive() : false;
  // If we haven't detected yet, nothing is incompatible
  if (autoDetect) { incompatPrior = salmon::math::LOG_1; }

  auto expectedLibraryFormat = readLib.format();
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

  auto isUnexpectedOrphan = [expectedLibraryFormat](AlnT& aln) -> bool {
    return (expectedLibraryFormat.type == ReadType::PAIRED_END and aln.mateStatus != rapmap::utils::MateStatus::PAIRED_END_PAIRED);
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

      //extract barcode of the read
      uint32_t barcode = alnGroup.barcode();
      uint64_t umi = alnGroup.umi();

      // We start out with probability 0
      double sumOfAlignProbs{LOG_0};

      // Record whether or not this read is unique to a single transcript.
      bool transcriptUnique{true};

      auto firstTranscriptID = alnGroup.alignments().front().transcriptID();
      std::unordered_set<size_t> observedTranscripts;

      std::vector<uint32_t> txpIDs;
      std::vector<double> auxProbs;
      double auxDenom = salmon::math::LOG_0;

      uint32_t numInGroup{0};
      uint32_t prevTxpID{0};

      hasCompatibleMapping = false;

      // For each alignment of this read
      for (auto& aln : alnGroup.alignments()) {

        useAuxParams = ((localNumAssignedFragments + numAssignedFragments) >= salmonOpts.numPreBurninFrags);
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
          if (isUnexpectedOrphan(aln)) {
            logFragProb = LOG_EPSILON;
          }

          if (flen > 0.0 and useFragLengthDist and considerCondProb) {
            size_t fl = flen;
            double lenProb = fragLengthDist.pmf(fl);
            if (burnedIn) {
              /* condition fragment length prob on txp length */
              double refLengthCM = fragLengthDist.cmf(static_cast<size_t>(refLength));
              bool computeMass = fl < refLength and !salmon::math::isLog0(refLengthCM);
              logFragProb = (computeMass) ?
                                      (lenProb - refLengthCM) :
                salmon::math::LOG_EPSILON;
              if (computeMass and refLengthCM < lenProb) {
                // Threading is hard!  It's possible that an update to the PMF
                //snuck in between when we asked to cache the CMF and when the
                // "burnedIn" variable was last seen as false.
                log->info("reference length = {}, CMF[refLen] = {}, fragLen = {},"
                          "PMF[fragLen] = {}", refLength, std::exp(refLengthCM),
                          aln.fragLength(), std::exp(lenProb));
              }
            } else if (useAuxParams) {
              logFragProb = lenProb;
            }
            //logFragProb = lenProb;
            //logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln.fragLength()));
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
          bool isCompat =
            salmon::utils::isCompatible(
                                        aln.libFormat(),
                                        expectedLibraryFormat,
                                        static_cast<int32_t>(aln.pos),
                                        aln.fwd,
                                        aln.mateStatus);
          double logAlignCompatProb = isCompat ? LOG_1 : incompatPrior;
          if (!isCompat and salmonOpts.ignoreIncompat) {
            aln.logProb = salmon::math::LOG_0;
            continue;
          }

          // Allow for a non-uniform fragment start position distribution

          double startPosProb{-logRefLength};
          // DEC 9
          if (aln.mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED and !noLengthCorrection) {
            startPosProb = (flen <= refLength) ? -std::log(refLength - flen + 1) : salmon::math::LOG_EPSILON;
          }

          double fragStartLogNumerator{salmon::math::LOG_1};
          double fragStartLogDenominator{salmon::math::LOG_1};

          auto hitPos = aln.hitPos();

          // Increment the count of this type of read that we've seen
          ++libTypeCounts[aln.libFormat().formatID()];
          //
          if (!hasCompatibleMapping and logAlignCompatProb == LOG_1) { hasCompatibleMapping = true; }

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

          if (updateCounts and
              observedTranscripts.find(transcriptID) ==
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
        if (hasCompatibleMapping) { ++numCompatibleFragments; }
      }

      auto eqSize = txpIDs.size();
      if (eqSize > 0) {
        if (useRankEqClasses and eqSize > 1) {
          std::vector<int> inds(eqSize);
          std::iota(inds.begin(), inds.end(), 0);
          // Get the indices in order by conditional probability
          std::sort(inds.begin(), inds.end(),
                    [&auxProbs](int i, int j) -> bool { return auxProbs[i] < auxProbs[j]; });
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

        TranscriptGroup tg(txpIDs);
        eqBuilder.addBarcodeGroup(std::move(tg), auxProbs, barcode, umi);
      }

      //+++++++++++++++++++++++++++++++++++++++
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

        double r = uni(randEng);
        if (!burnedIn and r < std::exp(aln.logProb)) {

          //Old fragment length calc: double fragLength = aln.fragLength();
          auto fragLength = aln.fragLengthPedantic(transcript.RefLength);
          if (fragLength > 0) {
            fragLengthDist.addVal(fragLength, logForgettingMass);
          }

          if (useFSPD) {
            auto hitPos = aln.hitPos();
            auto& fragStartDist = fragStartDists[transcript.lengthClassIndex()];
            fragStartDist.addVal(hitPos, transcript.RefLength,
                                 logForgettingMass);
          }
        }
      } // end normalize

      //+++++++++++++++++++++++++++++++++++++++

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
    } // end read group
  }   // end timer

  if (zeroProbFrags > 0) {
    auto batchReads = batchHits.size();
    maxZeroFrac = std::max(maxZeroFrac, static_cast<double>(100.0 * zeroProbFrags) / batchReads);
  }

  numAssignedFragments += localNumAssignedFragments;
  if (numAssignedFragments >= numBurninFrags and !burnedIn) {
    if (useFSPD) {
      // update all of the fragment start position
      // distributions
      for (auto& fspd : fragStartDists) {
        fspd.update();
      }
    }
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

/// START QUASI

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename RapMapIndexT, typename ProtocolT>
void processReadsQuasi(
                       paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
                       AlnGroupVec<SMEMAlignment>& structureVec,
                       std::atomic<uint64_t>& numObservedFragments,
                       std::atomic<uint64_t>& numAssignedFragments,
                       std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
                       RapMapIndexT* idx, std::vector<Transcript>& transcripts,
                       ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
                       FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
                       SalmonOpts& salmonOpts,
                       std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
                       volatile bool& writeToCache, AlevinOpts<ProtocolT>& alevinOpts,
                       SoftMapT& barcodeMap, spp::sparse_hash_map<std::string, uint32_t>& trBcs) {

  // ERROR
  salmonOpts.jointLog->error("MEM-mapping cannot be used with the Quasi index "
                             "--- please report this bug on GitHub");
  std::exit(1);
}

template <typename RapMapIndexT, typename ProtocolT>
void processReadsQuasi(
                       paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
                       AlnGroupVec<QuasiAlignment>& structureVec,
                       std::atomic<uint64_t>& numObservedFragments,
                       std::atomic<uint64_t>& numAssignedFragments,
                       std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
                       RapMapIndexT* qidx, std::vector<Transcript>& transcripts,
                       ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
                       FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
                       SalmonOpts& salmonOpts,
                       std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
                       volatile bool& writeToCache, AlevinOpts<ProtocolT>& alevinOpts,
                       SoftMapT& barcodeMap,
                       spp::sparse_hash_map<std::string, uint32_t>& trBcs) {
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
  spdlog::logger* unmappedLogger = (writeUnmapped) ? salmonOpts.unmappedLog.get() : nullptr;

  // Write unmapped reads
  fmt::MemoryWriter orphanLinks;
  bool writeOrphanLinks = salmonOpts.writeOrphanLinks;
  spdlog::logger* orphanLinkLogger = (writeOrphanLinks) ? salmonOpts.orphanLinkLog.get() : nullptr;

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
  bool strictIntersect = salmonOpts.strictIntersect;
  bool consistentHits = salmonOpts.consistentHits;
  bool quiet = salmonOpts.quiet;

  bool tooManyHits{false};
  size_t maxNumHits{salmonOpts.maxReadOccs};
  size_t readLenLeft{0};
  size_t readLenRight{0};
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

  SASearcher<RapMapIndexT> saSearcher(qidx);
  std::vector<QuasiAlignment> leftHits;
  std::vector<QuasiAlignment> rightHits;
  rapmap::utils::HitCounters hctr;
  salmon::utils::MappingType mapType{salmon::utils::MappingType::UNMAPPED};


  PairAlignmentFormatter<RapMapIndexT*> formatter(qidx);
  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);

  auto rg = parser->getReadGroup();
  while (parser->refill(rg)) {
    rangeSize = rg.size();

    if (rangeSize > structureVec.size()) {
      salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} "
                                 "--- this shouldn't happen.\n"
                                 "Please report this bug on GitHub",
                                 rangeSize, structureVec.size());
      std::exit(1);
    }

    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      auto& rp = rg[i];
      readLenLeft = rp.first.seq.length();
      readLenRight= rp.second.seq.length();

      bool tooShortLeft = (readLenLeft < minK);
      bool tooShortRight = (readLenRight < minK);
      tooManyHits = false;
      //localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();
      auto& jointHits = jointHitGroup.alignments();
      leftHits.clear();
      rightHits.clear();
      mapType = salmon::utils::MappingType::UNMAPPED;

      //////////////////////////////////////////////////////////////
      // extracting barcodes
      size_t barcodeLength = alevinOpts.protocol.barcodeLength;
      size_t umiLength = alevinOpts.protocol.umiLength;
      std::string umi, barcode;
      uint32_t barcodeIdx{0};
      bool lh, rh, isExtractOk, seqOk;

      if (alevinOpts.protocol.end == bcEnd::FIVE){
        if(alevinOpts.nobarcode){
          barcode = "AAA";
          isExtractOk = true;
          seqOk = true;
          alevinOpts.protocol.barcodeLength = 0;
        }
        else{
          isExtractOk = aut::extractBarcode(rp.first.seq,
                                            alevinOpts.protocol,
                                            barcode);
          seqOk = aut::sequenceCheck(barcode,
                                     alevinOpts,
                                     iomutex,
                                     Sequence::BARCODE);
        }

        //std::vector<std::string>::iterator bcIt = find(trBcs.begin(),
        //                                               trBcs.end(),
        //                                               barcode);

        //bool inTr = bcIt != trBcs.end();
        bool inTr = trBcs.contains(barcode);
        bool indOk {false};

        if(inTr){
          barcodeIdx = trBcs[barcode];
        }
        else{
          indOk = barcodeMap.find(barcode) != barcodeMap.end();
          if (indOk){
            barcode = barcodeMap[barcode].front().first;
            //auto trIt =  find(trBcs.begin(), trBcs.end(), barcode);
            bool trIt =  trBcs.contains(barcode);
            if(not trIt){
              salmonOpts.jointLog->error("Wrong entry in barcode softmap.\n"
                                         "Please Report this on github");
              exit(1);
            }
            else{
              barcodeIdx = trBcs[barcode];
            }
          }
        }

        if (not isExtractOk or (not inTr and not indOk) or not seqOk) {
          lh = rh = false;
        }
        else{
          //corrBarcodeIndex = barcodeMap[barcodeIndex];
          jointHitGroup.setBarcode(barcodeIdx);

          aut::extractUMI(rp.first.seq, alevinOpts.protocol, umi);
          alevin::kmer::AlvKmer umiIdx(umiLength);
          bool isUmiIdxOk = umiIdx.fromStr(umi);

          if(not isUmiIdxOk){
            salmonOpts.jointLog->error("umi indexing of jellyfish failing.\n"
                                       "Please report on github");
            exit(1);
          }

          jointHitGroup.setUMI(umiIdx.umiWord());
          //clearing barcode string to use as false mate
          barcode.clear();
          lh = tooShortLeft ? false : hitCollector(barcode,
                                                   leftHits, saSearcher,
                                                   MateStatus::PAIRED_END_LEFT,
                                                   consistentHits);
          rh = tooShortRight ? false : hitCollector(rp.second.seq,
                                                    rightHits, saSearcher,
                                                    MateStatus::PAIRED_END_RIGHT,
                                                    consistentHits);
        }
      }
      //else if (alevinOpts.barcodeEnd == THREE
      //  //have to handle 3 prime reverse complement
      //  lh = tooShortLeft ? false : hitCollector(rp.first.seq,
      //                                           leftHits, saSearcher,
      //                                           MateStatus::PAIRED_END_LEFT,
      //                                           consistentHits);

      //  rh = tooShortRight ? false : hitCollector(rp.second.seq,
      //                                            rightHits, saSearcher,
      //                                            MateStatus::PAIRED_END_RIGHT,
      //                                            consistentHits);
      //}
      else{
        salmonOpts.jointLog->error( "wrong barcode-end parameters.\n"
                                    "Please report this bug on Github");
        salmonOpts.jointLog->flush();
        spdlog::drop_all();
        std::exit(1);
      }
      //////////////////////////////////////////////////////////////
      // Consider a read as too short if both ends are too short
      if (tooShortLeft and tooShortRight) {
        ++shortFragStats.numTooShort;
        shortFragStats.shortest = std::min(shortFragStats.shortest,
                                           std::max(readLenLeft, readLenRight));
      } else {
        // If we actually attempted to map the fragment (it wasn't too short),
        // then
        // do the intersection.
        if (strictIntersect) {
          rapmap::utils::mergeLeftRightHits(leftHits, rightHits, jointHits,
                                            readLenLeft, maxNumHits,
                                            tooManyHits, hctr);
        } else {
          rapmap::utils::mergeLeftRightHitsFuzzy(lh, rh, leftHits, rightHits,
                                                 jointHits, readLenLeft,
                                                 maxNumHits, tooManyHits, hctr);
        }

        if (initialRound) {
          upperBoundHits += (jointHits.size() > 0);
        }

        // If the read mapped to > maxReadOccs places, discard it
        if (jointHits.size() > salmonOpts.maxReadOccs) {
          jointHitGroup.clearAlignments();
        }
      }

      // NOTE: This will currently not work with "strict intersect", i.e.
      // nothing will be output here with strict intersect.
      if (writeOrphanLinks) {
        // If we are not using strict intersection, then joint hits
        // can only be zero when either:
        // 1) there are *no* hits or
        // 2) there are hits for *both* the left and right reads, but not to the same txp
        if (!strictIntersect and jointHits.size() == 0) {
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
      bool isPaired{false};
      if (jointHits.size() > 0) {
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
        } else {
          // If these aren't paired-end reads --- so that
          // we have orphans --- make sure we sort the
          // mappings so that they are in transcript order
          if (!isPaired) {
            // Find the end of the hits for the left read
            auto leftHitEndIt = std::partition_point(
                                                     jointHits.begin(), jointHits.end(),
                                                     [](const QuasiAlignment& q) -> bool {
                                                       return q.mateStatus ==
                                                         rapmap::utils::MateStatus::PAIRED_END_LEFT;
                                                     });
            // If we found left hits
            bool foundLeftMappings = (leftHitEndIt > jointHits.begin());
            // If we found right hits
            bool foundRightMappings = (leftHitEndIt  < jointHits.end());

            if (foundLeftMappings and foundRightMappings) {
              mapType = salmon::utils::MappingType::BOTH_ORPHAN;
            } else if (foundLeftMappings) {
              mapType = salmon::utils::MappingType::LEFT_ORPHAN;
            } else if (foundRightMappings) {
              mapType = salmon::utils::MappingType::RIGHT_ORPHAN;
            }

            // Merge the hits so that the entire list is in order
            // by transcript ID.
            std::inplace_merge(
                               jointHits.begin(), leftHitEndIt, jointHits.end(),
                               [](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
                                 return a.transcriptID() < b.transcriptID();
                               });
          }
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
                   startPos2 + readBias2.contextAfter(read2RC) <
                   static_cast<int32_t>(t.RefLength))) {

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
          case MateStatus::SINGLE_END : {
            // do nothing
            //salmonOpts.jointLog->warn("");
          } break;
          }
        }

        if (writeQuasimappings) {
          rapmap::utils::writeAlignmentsToStream(rp, formatter,
                                                 hctr, jointHits, sstream);
        }

      } else {
        // This read was completely unmapped.
        mapType = salmon::utils::MappingType::UNMAPPED;
      }

      if (writeUnmapped and mapType != salmon::utils::MappingType::PAIRED_MAPPED) {
        // If we have no mappings --- then there's nothing to do
        // unless we're outputting names for un-mapped reads
        unmappedNames << rp.first.name << ' ' << salmon::utils::str(mapType) << '\n';
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
          fmt::print(stderr, "\033[A\r\r{}processed{} {} Million {}fragments{}\n",
                     green, red, numObservedFragments/1000000, green, RESET_COLOR);
          fmt::print(stderr, "hits: {}, hits per frag:  {}", validHits,
                     validHits / static_cast<float>(prevObservedFrags));
        } else {
          fmt::print(stderr, "\r\r{}processed{} {} {}fragments{}", green, red,
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
    AlnGroupVecRange<QuasiAlignment> hitLists = boost::make_iterator_range(
                                                                           structureVec.begin(), structureVec.begin() + rangeSize);
    processMiniBatch<QuasiAlignment>(
                                     readExp, fmCalc, firstTimestepOfRound, rl, salmonOpts, hitLists,
                                     transcripts, clusterForest, fragLengthDist, observedBiasParams,
                                     numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac,
                                     alevinOpts, iomutex);
  }

  if (maxZeroFrac > 0.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }

  readExp.updateShortFrags(shortFragStats);
}

template <typename AlnT, typename ProtocolT>
void processReadLibrary(
                        ReadExperimentT& readExp, ReadLibrary& rl, SalmonIndex* sidx,
                        std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
                        std::atomic<uint64_t>&
                        numObservedFragments, // total number of reads we've looked at
                        std::atomic<uint64_t>&
                        numAssignedFragments,              // total number of assigned reads
                        std::atomic<uint64_t>& upperBoundHits, // upper bound on # of mapped frags
                        bool initialRound,
                        std::atomic<bool>& burnedIn, ForgettingMassCalculator& fmCalc,
                        FragmentLengthDistribution& fragLengthDist,
                        SalmonOpts& salmonOpts, bool greedyChain,
                        std::mutex& iomutex, size_t numThreads,
                        std::vector<AlnGroupVec<AlnT>>& structureVec, volatile bool& writeToCache,
                        AlevinOpts<ProtocolT>& alevinOpts,
                        SoftMapT& barcodeMap,
                        spp::sparse_hash_map<std::string, uint32_t>& trBcs) {

  std::vector<std::thread> threads;

  std::atomic<uint64_t> numValidHits{0};
  rl.checkValid();

  auto indexType = sidx->indexType();

  std::unique_ptr<paired_parser> pairedParserPtr{nullptr};

  /** sequence-specific and GC-fragment bias vectors --- each thread gets it's
   * own **/
  std::vector<BiasParams> observedBiasParams(numThreads,
                                             BiasParams(salmonOpts.numConditionalGCBins, salmonOpts.numFragGCBins, false));

  // If the read library is paired-end
  // ------ Paired-end --------
  if (rl.format().type == ReadType::PAIRED_END) {

    if (rl.mates1().size() != rl.mates2().size()) {
      salmonOpts.jointLog->error("The number of provided files for "
                                 "-1 and -2 must be the same!");
      std::exit(1);
    }

    size_t numFiles = rl.mates1().size() + rl.mates2().size();
    uint32_t numParsingThreads{1};
    // HACK!
    if(numThreads > 1){
      numThreads -= 1;
    }
    if (rl.mates1().size() > 1 and numThreads > 8) { numParsingThreads = 2; numThreads -= 1;}
    pairedParserPtr.reset(new paired_parser(rl.mates1(), rl.mates2(), numThreads, numParsingThreads, miniBatchSize));
    pairedParserPtr->start();

    // True if we have a 64-bit SA index, false otherwise
    bool largeIndex = sidx->is64BitQuasi();
    bool perfectHashIndex = sidx->isPerfectHashQuasi();
    for (size_t i = 0; i < numThreads; ++i) {
      // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
      // change value before the lambda below is evaluated --- crazy!
      if (largeIndex) {
        if (perfectHashIndex) { // Perfect Hash
          if (salmonOpts.qmFileName != "" and i == 0) {
            rapmap::utils::writeSAMHeader(*(sidx->quasiIndexPerfectHash64()), salmonOpts.qmLog);
          }
          auto threadFun = [&, i]() -> void {
            processReadsQuasi<RapMapSAIndex<int64_t, PerfectHash<int64_t>>>(
                                                                            pairedParserPtr.get(), readExp, rl, structureVec[i],
                                                                            numObservedFragments, numAssignedFragments, numValidHits,
                                                                            upperBoundHits, sidx->quasiIndexPerfectHash64(), transcripts,
                                                                            fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                                                                            salmonOpts, iomutex, initialRound,
                                                                            burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs);
          };
          threads.emplace_back(threadFun);
        } else { // Dense Hash
          if (salmonOpts.qmFileName != "" and i == 0) {
            rapmap::utils::writeSAMHeader(*(sidx->quasiIndex64()), salmonOpts.qmLog);
          }
          auto threadFun = [&, i]() -> void {
            processReadsQuasi<RapMapSAIndex<int64_t, DenseHash<int64_t>>>(
                                                                          pairedParserPtr.get(), readExp, rl, structureVec[i],
                                                                          numObservedFragments, numAssignedFragments, numValidHits,
                                                                          upperBoundHits, sidx->quasiIndex64(), transcripts, fmCalc,
                                                                          clusterForest, fragLengthDist, observedBiasParams[i],
                                                                          salmonOpts, iomutex, initialRound,
                                                                          burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs);
          };
          threads.emplace_back(threadFun);
        }
      } else {
        if (perfectHashIndex) { // Perfect Hash
          if (salmonOpts.qmFileName != "" and i == 0) {
            rapmap::utils::writeSAMHeader(*(sidx->quasiIndexPerfectHash32()), salmonOpts.qmLog);
          }
          auto threadFun = [&, i]() -> void {
            processReadsQuasi<RapMapSAIndex<int32_t, PerfectHash<int32_t>>>(
                                                                            pairedParserPtr.get(), readExp, rl, structureVec[i],
                                                                            numObservedFragments, numAssignedFragments, numValidHits,
                                                                            upperBoundHits, sidx->quasiIndexPerfectHash32(), transcripts,
                                                                            fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                                                                            salmonOpts, iomutex, initialRound,
                                                                            burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs);
          };
          threads.emplace_back(threadFun);
        } else { // Dense Hash
          if (salmonOpts.qmFileName != "" and i == 0) {
            rapmap::utils::writeSAMHeader(*(sidx->quasiIndex32()), salmonOpts.qmLog);
          }
          auto threadFun = [&, i]() -> void {
            processReadsQuasi<RapMapSAIndex<int32_t, DenseHash<int32_t>>>(
                                                                          pairedParserPtr.get(), readExp, rl, structureVec[i],
                                                                          numObservedFragments, numAssignedFragments, numValidHits,
                                                                          upperBoundHits, sidx->quasiIndex32(), transcripts, fmCalc,
                                                                          clusterForest, fragLengthDist, observedBiasParams[i],
                                                                          salmonOpts, iomutex, initialRound,
                                                                          burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs);
          };
          threads.emplace_back(threadFun);
        }

      } // End spawn current thread

    } // End spawn all threads

    for (auto& t : threads) {
      t.join();
    }

    pairedParserPtr->stop();

    //+++++++++++++++++++++++++++++++++++++++
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

      auto& fw = readExp.readBiasModelObserved(salmon::utils::Direction::FORWARD);
      auto& rc =
        readExp.readBiasModelObserved(salmon::utils::Direction::REVERSE_COMPLEMENT);

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

    //+++++++++++++++++++++++++++++++++++++++

  }
}

/**
 *  Quantify the targets given in the file `transcriptFile` using the
 *  reads in the given set of `readLibraries`, and write the results
 *  to the file `outputFile`.  The reads are assumed to be in the format
 *  specified by `libFmt`.
 *
 */
template <typename AlnT, typename ProtocolT>
void quantifyLibrary(ReadExperimentT& experiment, bool greedyChain,
                     SalmonOpts& salmonOpts,
                     uint32_t numQuantThreads,
                     AlevinOpts<ProtocolT>& alevinOpts,
                     SoftMapT& barcodeMap,
                     spp::sparse_hash_map<std::string, uint32_t>& trBcs) {

  bool burnedIn = (salmonOpts.numBurninFrags == 0);
  uint64_t numRequiredFragments = salmonOpts.numRequiredFragments;
  std::atomic<uint64_t> upperBoundHits{0};
  // ErrorModel errMod(1.00);
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
                             greedyChain, ioMutex,
                             numQuantThreads, groupVec, writeToCache,
                             alevinOpts, barcodeMap, trBcs);

    numAssignedFragments = totalAssignedFragments - prevNumAssignedFragments;
  };

  if (!salmonOpts.quiet) {
    salmonOpts.jointLog->flush();
    fmt::print(stderr, "\n\n\n\n");
  }

  // Process all of the reads
  experiment.processReads(numQuantThreads, salmonOpts,
                          processReadLibraryCallback);
  experiment.setNumObservedFragments(numObservedFragments);

  // EQCLASS
  // changing it to alevin based finish
  bool done = experiment.equivalenceClassBuilder().alv_finish();
  // skip the extra online rounds

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

  //+++++++++++++++++++++++++++++++++++++++
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

    // If we didn't have a sufficient number of samples for burnin,
    // then also ignore modeling of the fragment start position
    // distribution.
    if (salmonOpts.useFSPD) {
      salmonOpts.useFSPD = false;
      jointLog->warn("Since only {} (< {}) fragments were observed, modeling "
                     "of the fragment start position "
                     "distribution has been disabled",
                     totalAssignedFragments, salmonOpts.numBurninFrags);
    }
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
    if (salmonOpts.allowOrphans) {
      double mappingRate = totalAssignedFragments.load() /
        static_cast<double>(numObservedFragments.load());
      experiment.setEffectiveMappingRate(mappingRate);
    } else {
      experiment.setEffectiveMappingRate(upperBoundMappingRate);
    }
  }
  //+++++++++++++++++++++++++++++++++++++++
  jointLog->info("Mapping rate = {}\%\n",
                 experiment.effectiveMappingRate() * 100.0);
  jointLog->info("finished quantifyLibrary()");
}


template <typename ProtocolT>
int alevinQuant(AlevinOpts<ProtocolT>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter, size_t numLowConfidentBarcode){

  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;
  try{
    //auto fileLog = sopt.fileLog;
    auto jointLog = aopt.jointLog;
    auto indexDirectory = sopt.indexDirectory;
    auto outputDirectory = sopt.outputDirectory;
    bool greedyChain = true;

    jointLog->info("parsing read library format");

    // ==== Library format processing ===
    vector<ReadLibrary> readLibraries =
      salmon::utils::extractReadLibraries(orderedOptions);

    if (readLibraries.size() == 0) {
      jointLog->error("Failed to successfully parse any complete read libraries."
                      " Please make sure you provided arguments properly to -1, -2 (for paired-end libraries)"
                      " or -r (for single-end libraries), and that the library format option (-l) comes before,"
                      " the read libraries.");
      std::exit(1);
    }
    // ==== END: Library format processing ===

    SalmonIndexVersionInfo versionInfo;
    boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
    versionInfo.load(versionPath);
    auto idxType = versionInfo.indexType();

    ReadExperimentT experiment(readLibraries, indexDirectory, sopt);

    // This will be the class in charge of maintaining our
    // rich equivalence classes
    experiment.equivalenceClassBuilder().start();

    auto indexType = experiment.getIndex()->indexType();

    // We can only do fragment GC bias correction, for the time being, with
    // paired-end reads
    if (sopt.gcBiasCorrect) {
      for (auto& rl : readLibraries) {
        if (rl.format().type != ReadType::PAIRED_END) {
          jointLog->warn("Fragment GC bias correction is currently *experimental* "
                         "in single-end libraries.  Please use this option "
                         "with caution.");
          //sopt.gcBiasCorrect = false;
        }
      }
    }

    std::vector<std::string> trueBarcodesVec (trueBarcodes.begin(), trueBarcodes.end());
    spp::sparse_hash_map<std::string, uint32_t> trueBarcodesIndexMap;
    for(size_t i=0; i<trueBarcodes.size(); i++){
      trueBarcodesIndexMap[ trueBarcodesVec[i] ] = i;
    }

    sopt.allowOrphans = true;
    sopt.useQuasi = true;
    if(sopt.numThreads > 1){
      sopt.numThreads -= 1;
    }
    quantifyLibrary<QuasiAlignment>(experiment, greedyChain, sopt,
                                    sopt.numThreads, aopt,
                                    barcodeMap, trueBarcodesIndexMap);

    // Write out information about the command / run
    salmon::utils::writeCmdInfo(sopt, orderedOptions);

    GZipWriter gzw(outputDirectory, jointLog);
    //+++++++++++++++++++++++++++++++++++++++
    // Now that the streaming pass is complete, we have
    // our initial estimates, and our rich equivalence
    // classes.  Perform further optimization until
    // convergence.
    // NOTE: A side-effect of calling the optimizer is that
    // the `EffectiveLength` field of each transcript is
    // set to its final value.
    if(aopt.dumpBarcodeEq){
      gzw.writeEquivCounts(aopt, experiment);
    }

    std::vector<uint32_t> umiCount(trueBarcodesVec.size());
    auto& fullEqMap = experiment.equivalenceClassBuilder().eqMap();
    for(auto& eq: fullEqMap.lock_table()){
      auto& bg = eq.second.barcodeGroup;
      for(auto& bcIt: bg){
        size_t bcCount{0};
        for(auto& ugIt: bcIt.second){
          bcCount += ugIt.second;
        }
        auto bc = bcIt.first;
        umiCount[bc] += bcCount;
      }
    }

    if(aopt.dumpfeatures){
      auto mapCountFileName = aopt.outputDirectory / "MappedUmi.txt";
      std::ofstream mFile;
      mFile.open(mapCountFileName.string());

      for(size_t i=0; i<umiCount.size(); i++){
        mFile<<trueBarcodesVec[i]<< "\t"<<umiCount[i]<<"\n";
      }
      mFile.close();
    }
    ////////////////////////////////////////////
    // deduplication starts from here
    ////////////////////////////////////////////

    if(aopt.dedup) {
      jointLog->info("Starting optimizer");
      jointLog->flush();
      std::cout<< "\n\n";

      CollapsedCellOptimizer optimizer;
      bool optSuccess = optimizer.optimize(experiment, aopt,
                                           gzw,
                                           trueBarcodesVec,
                                           umiCount,
                                           freqCounter,
                                           numLowConfidentBarcode);
      if (!optSuccess) {
        jointLog->error(
                        "The optimization algorithm failed. This is likely the result of "
                        "bad input (or a bug). If you cannot track down the cause, please "
                        "report this issue on GitHub.");
        jointLog->flush();
        exit(1);
      }
      std::cout<< "\n\n";
      jointLog->info("Finished optimizer");
    }
    else{
      jointLog->warn("No Dedup command given, is it what you want?");
    }

    jointLog->flush();

    //CollapsedJointEMOptimizer optimizer;
    //salmon::utils::normalizeAlphas(sopt, experiment);

    //TrueMapT trBcMap;
    //uint32_t ind{0};
    //for(const auto& trBc : trueBarcodes){
    //  trBcMap[trBc] = ind;
    //  ind += 1;
    //}
    //bool optSuccess = optimizer.optimize(experiment, sopt,
    //                                     barcodeMap, trBcMap,
    //                                     0.01, 10000);

    //if (!optSuccess) {
    //  jointLog->error(
    //                  "The optimization algorithm failed. This is likely the result of "
    //                  "bad input (or a bug). If you cannot track down the cause, please "
    //                  "report this issue on GitHub.");
    //  return 1;
    //}
    //jointLog->info("Finished optimizer");

    //free(memOptions);
    //size_t tnum{0};

    //jointLog->info("writing output \n");

    //bfs::path estFilePath = outputDirectory / "quant.sf";

    // Write the main results
    //gzw.writeAbundances(sopt, experiment);

    // If we are dumping the equivalence classes, then
    // do it here.
    //if (sopt.dumpEq) {
    //  gzw.writeEquivCounts(sopt, experiment);
    //}

    //if (sopt.numGibbsSamples > 0) {

    //  jointLog->info("Starting Gibbs Sampler");
    //  CollapsedGibbsSampler sampler;
    //  gzw.setSamplingPath(sopt);
    //  // The function we'll use as a callback to write samples
    //  std::function<bool(const std::vector<double>&)> bsWriter =
    //    [&gzw](const std::vector<double>& alphas) -> bool {
    //    return gzw.writeBootstrap(alphas, true);
    //  };

    //  bool sampleSuccess =
    //    //sampler.sampleMultipleChains(experiment, sopt, bsWriter, sopt.numGibbsSamples);
    //    sampler.sample(experiment, sopt, bsWriter, sopt.numGibbsSamples);
    //  if (!sampleSuccess) {
    //    jointLog->error("Encountered error during Gibbs sampling.\n"
    //                    "This should not happen.\n"
    //                    "Please file a bug report on GitHub.\n");
    //    return 1;
    //  }
    //  jointLog->info("Finished Gibbs Sampler");
    //} else if (sopt.numBootstraps > 0) {
    //  gzw.setSamplingPath(sopt);
    //  // The function we'll use as a callback to write samples
    //  std::function<bool(const std::vector<double>&)> bsWriter =
    //    [&gzw](const std::vector<double>& alphas) -> bool {
    //    return gzw.writeBootstrap(alphas);
    //  };

    //  jointLog->info("Starting Bootstrapping");
    //  bool bootstrapSuccess =
    //    optimizer.gatherBootstraps(experiment, sopt, bsWriter, 0.01, 10000);
    //  jointLog->info("Finished Bootstrapping");
    //  if (!bootstrapSuccess) {
    //    jointLog->error("Encountered error during bootstrapping.\n"
    //                    "This should not happen.\n"
    //                    "Please file a bug report on GitHub.\n");
    //    return 1;
    //  }
    //}
    //+++++++++++++++++++++++++++++++++++++++

    bfs::path libCountFilePath = outputDirectory / "lib_format_counts.json";
    experiment.summarizeLibraryTypeCounts(libCountFilePath);

    //+++++++++++++++++++++++++++++++++++++++
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

    /** If the user requested gene-level abundances, then compute those now **/
    //if (vm.count("geneMap")) {
    //  try {
    //    salmon::utils::generateGeneLevelEstimates(sopt.geneMapPath,
    //                                              outputDirectory);
    //  } catch (std::invalid_argument& e) {
    //    fmt::print(stderr, "Error: [{}] when trying to compute gene-level "
    //               "estimates. The gene-level file(s) may not exist",
    //               e.what());
    //  }
    //}

    //+++++++++++++++++++++++++++++++++++++++


    if (sopt.writeUnmappedNames) {
      auto l = sopt.unmappedLog.get();
      // If the logger was created, then flush it and
      // close the associated file.
      if (l) {
        l->flush();
        if (sopt.unmappedFile) { sopt.unmappedFile->close(); }
      }
    }

    if (sopt.writeOrphanLinks) {
      auto l = sopt.orphanLinkLog.get();
      // If the logger was created, then flush it and
      // close the associated file.
      if (l) {
        l->flush();
        if (sopt.orphanLinkFile) { sopt.orphanLinkFile->close(); }
      }
    }

    // if we wrote quasimappings, flush that buffer
    if (sopt.qmFileName != "" ){
      sopt.qmLog->flush();
      // if we wrote to a buffer other than stdout, close
      // the file
      if (sopt.qmFileName != "-") { sopt.qmFile.close(); }
    }

    sopt.runStopTime = salmon::utils::getCurrentTimeAsString();

    // Write meta-information about the run
    gzw.writeMeta(sopt, experiment);

  } catch (po::error& e) {
    std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << " alevin was invoked improperly.\n";
    std::cerr << "For usage information, try "
              << " alevin --help\nExiting.\n";
    std::exit(1);
  }

  return 0;
}

namespace apt = alevin::protocols;
template
int alevinQuant(AlevinOpts<apt::DropSeq>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::InDrop>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::Chromium>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::Custom>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);

