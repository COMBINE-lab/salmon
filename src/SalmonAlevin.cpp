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

// core includes
#include "core/range.hpp"

//alevin include
#include "AlevinOpts.hpp"
#include "AlevinUtils.hpp"
#include "SingleCellProtocols.hpp"
#include "CollapsedCellOptimizer.hpp"
#include "BarcodeGroup.hpp"

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
#include "BiasParams.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "ForgettingMassCalculator.hpp"
#include "FragmentLengthDistribution.hpp"
#include "GZipWriter.hpp"
#include "HitManager.hpp"

#include "RapMapUtils.hpp"
#include "ReadExperiment.hpp"
#include "SACollector.hpp"
#include "SASearcher.hpp"
#include "SalmonOpts.hpp"
#include "PairAlignmentFormatter.hpp"
#include "SingleAlignmentFormatter.hpp"
#include "RapMapUtils.hpp"
#include "ksw2pp/KSW2Aligner.hpp"
#include "tsl/hopscotch_map.h"
#include "SelectiveAlignmentUtils.hpp"

namespace alevin{

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

  using CellBarcodeT = uint32_t;
  using UMIBarcodeT = uint64_t;

  template <typename AlnT> using AlevinAlnGroup = AlignmentGroup<AlnT, CellBarcodeT, UMIBarcodeT>;
  template <typename AlnT> using AlnGroupVec = std::vector<AlevinAlnGroup<AlnT>>;

  template <typename AlnT>
  using AlnGroupVecRange = core::range<typename AlnGroupVec<AlnT>::iterator>;

#define __MOODYCAMEL__
#if defined(__MOODYCAMEL__)
  template <typename AlnT>
  using AlnGroupQueue = moodycamel::ConcurrentQueue<AlevinAlnGroup<AlnT>*>;
#else
  template <typename AlnT>
  using AlnGroupQueue = tbb::concurrent_queue<AlevinAlnGroup<AlnT>*>;
#endif

  //#include "LightweightAlignmentDefs.hpp"
}

//have to create new namespace because of multiple definition
using namespace alevin;

/* ALEVIN DECLERATIONS*/
using bcEnd = BarcodeEnd;
namespace aut = alevin::utils;
using BlockedIndexRange = tbb::blocked_range<size_t>;
using ReadExperimentT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;

//Note: redundant code starts, first used in SalmonQuantify.cpp
using AlnCacheMap = selective_alignment::utils::AlnCacheMap;
/////// REDUNDANT CODE END//

template <typename AlnT, typename ProtocolT>
void processMiniBatchSimple(ReadExperimentT& readExp, ForgettingMassCalculator& fmCalc,
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

  bool updateCounts = initialRound;
  double incompatPrior = salmonOpts.incompatPrior;
  bool useReadCompat = incompatPrior != salmon::math::LOG_1;

  // If we're auto detecting the library type
  auto* detector = readLib.getDetector();
  bool autoDetect = (detector != nullptr) ? detector->isActive() : false;
  // If we haven't detected yet, nothing is incompatible
  if (autoDetect) { incompatPrior = salmon::math::LOG_1; }

  double logForgettingMass{0.0};
  uint64_t currentMinibatchTimestep{0};

  // logForgettingMass and currentMinibatchTimestep are OUT parameters!
  fmCalc.getLogMassAndTimestep(logForgettingMass, currentMinibatchTimestep);

  auto expectedLibraryFormat = readLib.format();
  uint64_t zeroProbFrags{0};

  // EQClass
  auto& eqBuilder = readExp.equivalenceClassBuilder();

  // Build reverse map from transcriptID => hit id
  using HitID = uint32_t;

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
      std::vector<uint32_t> txpIDs;

      uint32_t numInGroup{0};
      uint32_t prevTxpID{0};

      hasCompatibleMapping = false;

      // For each alignment of this read
      for (auto& aln : alnGroup.alignments()) {
        auto transcriptID = aln.transcriptID();
        auto& transcript = transcripts[transcriptID];
        transcriptUnique =
          transcriptUnique and (transcriptID == firstTranscriptID);

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

          // The probability that the fragments align to the given strands in
          // the given orientations.
          bool isCompat =
            salmon::utils::isCompatible(
                                        aln.libFormat(),
                                        expectedLibraryFormat,
                                        static_cast<int32_t>(aln.pos),
                                        aln.fwd,
                                        aln.mateStatus);
          double logAlignCompatProb = isCompat ? LOG_1 : incompatPrior;
          aln.logProb = logAlignCompatProb;
          if (!isCompat and salmonOpts.ignoreIncompat) {
            aln.logProb = salmon::math::LOG_0;
            continue;
          }

          // Increment the count of this type of read that we've seen
          ++libTypeCounts[aln.libFormat().formatID()];
          if (!hasCompatibleMapping and logAlignCompatProb == LOG_1) { hasCompatibleMapping = true; }

          // If this alignment had a zero probability, then skip it
          if (std::abs(aln.logProb) == LOG_0) {
            continue;
          }

          sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);

          if (transcriptID < prevTxpID) {
            std::cerr << "[ERROR] Transcript IDs are not in sorted order; "
              "please report this bug on GitHub!\n";
          }
          prevTxpID = transcriptID;
          txpIDs.push_back(transcriptID);
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
        TranscriptGroup tg(txpIDs);
        eqBuilder.addBarcodeGroup(std::move(tg), barcode, umi);
      }
            // update the single target transcript
      if (transcriptUnique) {
        if (updateCounts) {
          transcripts[firstTranscriptID].addUniqueCount(1);
        }
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
template <typename RapMapIndexT, typename ProtocolT>
void processReadsQuasi(
                       paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
                       AlnGroupVec<QuasiAlignment>& structureVec,
                       std::atomic<uint64_t>& numObservedFragments,
                       std::atomic<uint64_t>& numAssignedFragments,
                       std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
                       std::atomic<uint32_t>& smallSeqs,
                       std::atomic<uint32_t>& nSeqs,
                       RapMapIndexT* qidx, std::vector<Transcript>& transcripts,
                       ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
                       FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
                       SalmonOpts& salmonOpts,
                       std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
                       volatile bool& writeToCache, AlevinOpts<ProtocolT>& alevinOpts,
                       SoftMapT& barcodeMap,
                       spp::sparse_hash_map<std::string, uint32_t>& trBcs,
                       MappingStatistics& mstats
                       /*,std::vector<uint64_t>& uniqueFLD*/) {
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
  bool consistentHits = salmonOpts.consistentHits;
  bool quiet = salmonOpts.quiet;

  size_t maxNumHits{salmonOpts.maxReadOccs};
  size_t readLenLeft{0};
  size_t readLenRight{0};
  SACollector<RapMapIndexT> hitCollector(qidx);

  rapmap::utils::MappingConfig mc;
  mc.consistentHits = consistentHits;
  mc.doChaining = salmonOpts.validateMappings;
  mc.consensusFraction = (salmonOpts.consensusSlack == 0.0) ? 1.0 : (1.0 - salmonOpts.consensusSlack);

  rapmap::hit_manager::HitCollectorInfo<rapmap::utils::SAIntervalHit<typename RapMapIndexT::IndexType>> hcInfo;

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

  SASearcher<RapMapIndexT> saSearcher(qidx);
  rapmap::utils::HitCounters hctr;
  salmon::utils::MappingType mapType{salmon::utils::MappingType::UNMAPPED};


  PairAlignmentFormatter<RapMapIndexT*> formatter(qidx);
  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);

  //////////////////////
  // NOTE: validation mapping based new parameters
  std::string rc1; rc1.reserve(300);

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
  int8_t a = salmonOpts.matchScore;
  int8_t b = salmonOpts.mismatchPenalty;
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

  size_t numDropped{0};
  size_t numDecoyFrags{0};
  std::string readSubSeq;

  //std::vector<salmon::mapping::CacheEntry> alnCache; alnCache.reserve(15);
  AlnCacheMap alnCache; alnCache.reserve(16);
  //////////////////////

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

      bool tooShortRight = (readLenRight < (minK+alevinOpts.trimRight));
      //localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();
      auto& jointHits = jointHitGroup.alignments();
      hcInfo.clear();
      mapType = salmon::utils::MappingType::UNMAPPED;

      //////////////////////////////////////////////////////////////
      // extracting barcodes
      size_t barcodeLength = alevinOpts.protocol.barcodeLength;
      size_t umiLength = alevinOpts.protocol.umiLength;
      std::string umi;//, barcode;
      nonstd::optional<std::string> barcode;
      nonstd::optional<uint32_t> barcodeIdx;
      bool seqOk;

      if (alevinOpts.protocol.end == bcEnd::FIVE){
        barcode = aut::extractBarcode(rp.first.seq, alevinOpts.protocol);
        seqOk = (barcode.has_value()) ?
          aut::sequenceCheck(*barcode, Sequence::BARCODE) : false;

        // If we have a barcode sequence, but not yet an index
        if (seqOk and (not barcodeIdx)) {
          // If we get here, we have a sequence-valid barcode.
          // Check if it is in the trBcs map.
          auto trIt = trBcs.find(*barcode);

          // If it is, use that index
          if(trIt != trBcs.end()){
            barcodeIdx = trIt->second;
          } else{
            // If it's not, see if it's in the barcode map
            auto indIt = barcodeMap.find(*barcode);
            // If so grab the representative and get its index
            if (indIt != barcodeMap.end()){
              barcode = indIt->second.front().first;
              auto trItLoc = trBcs.find(*barcode);
              if(trItLoc == trBcs.end()){
                salmonOpts.jointLog->error("Wrong entry in barcode softmap.\n"
                                           "Please Report this on github");
                salmonOpts.jointLog->flush();
                spdlog::drop_all();
                exit(1);
              } else{
                barcodeIdx = trItLoc->second;
              }
            }
            // If it wasn't in the barcode map, it's not valid
            // and we should leave barcodeIdx as nullopt.
          }
        }

        // If we have a valid barcode
        if (barcodeIdx) {
          //corrBarcodeIndex = barcodeMap[barcodeIndex];
          jointHitGroup.setBarcode(*barcodeIdx);
          aut::extractUMI(rp.first.seq, alevinOpts.protocol, umi);

          if ( umiLength != umi.size() ) {
            smallSeqs += 1;
          } else{
            alevin::types::AlevinUMIKmer umiIdx;
            bool isUmiIdxOk = umiIdx.fromChars(umi);

            if(isUmiIdxOk){
              jointHitGroup.setUMI(umiIdx.word(0));

              auto seq_len = rp.second.seq.size();
              if (alevinOpts.trimRight > 0) {
                if ( !tooShortRight ) {
                  //std::string sub_seq = rp.second.seq.substr(0, seq_len-alevinOpts.trimRight);
                  //auto rh = hitCollector(sub_seq, saSearcher, hcInfo);
                  readSubSeq = rp.second.seq.substr(0, seq_len-alevinOpts.trimRight);
                  auto rh = hitCollector(readSubSeq, saSearcher, hcInfo);
                }
              } else {
                readSubSeq = rp.second.seq;
                auto rh = tooShortRight ? false : hitCollector(rp.second.seq, saSearcher, hcInfo);
              }
              rapmap::hit_manager::hitsToMappingsSimple(*qidx, mc,
                                                        MateStatus::PAIRED_END_RIGHT,
                                                        hcInfo, jointHits);
            } else{
              nSeqs += 1;
            }
          }
        }
      } else{
        salmonOpts.jointLog->error( "wrong barcode-end parameters.\n"
                                    "Please report this bug on Github");
        salmonOpts.jointLog->flush();
        spdlog::drop_all();
        std::exit(1);
      }
      //////////////////////////////////////////////////////////////
      // Consider a read as too short if the ``non-barcode'' end is too short
      if (tooShortRight) {
        ++shortFragStats.numTooShort;
        shortFragStats.shortest = std::min(shortFragStats.shortest,
                                           std::max(readLenLeft, readLenRight));
      }

      if (initialRound) {
        upperBoundHits += (jointHits.size() > 0);
      }

      // If the read mapped to > maxReadOccs places, discard it
      if (jointHits.size() > salmonOpts.maxReadOccs) {
        jointHitGroup.clearAlignments();
      }

        // adding validate mapping code
        bool tryAlign{salmonOpts.validateMappings};
        if (tryAlign and !jointHits.empty()) {
          alnCache.clear();
          /*
          auto seq_len = rp.second.seq.size();
          std::string sub_seq = rp.second.seq.substr(0, seq_len-alevinOpts.trimRight);
          auto* r1 = sub_seq.data();
          auto l1 = static_cast<int32_t>(sub_seq.length());
          */
          auto* r1 = readSubSeq.data();
          auto l1 = static_cast<int32_t>(readSubSeq.length());

          char* r1rc = nullptr;
          int32_t bestScore{std::numeric_limits<int32_t>::min()};
          int32_t bestDecoyScore{std::numeric_limits<int32_t>::lowest()};
          std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
          std::vector<bool> decoyVec(jointHits.size(), false);
          size_t idx{0};
          double optFrac{salmonOpts.minScoreFraction};
          int32_t maxReadScore{a * static_cast<int32_t>(readSubSeq.length())};

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
              rapmap::utils::reverseRead(readSubSeq, rc1);
              // we will not break the const promise
              r1rc = const_cast<char*>(rc1.data());
            }

            auto* rptr = h.fwd ? r1 : r1rc;
            int32_t s =
              selective_alignment::utils::getAlnScore(aligner, ez, h.pos, rptr, l1, tseq, tlen, a, b,
                                                      maxReadScore, h.chainStatus.getLeft(),
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
                                           [&ctr, &scores, &decoyVec, &numDropped, bestScore] (const QuasiAlignment& qa) -> bool {
                                             // soft filter
                                             //bool rem = (scores[ctr] == std::numeric_limits<int32_t>::min());
                                             //strict filter
                                             bool rem = decoyVec[ctr] ? true : (scores[ctr] < bestScore);
                                             ++ctr;
                                             numDropped += rem ? 1 : 0;
                                             return rem;
                                           }),
                            jointHits.end()
                            );
            // soft filter
            double bestScoreD = static_cast<double>(bestScore);
            std::for_each(jointHits.begin(), jointHits.end(),
                          [bestScoreD, writeQuasimappings](QuasiAlignment& qa) -> void {
                            if (writeQuasimappings) { qa.alnScore(static_cast<int32_t>(qa.score())); }
                            double v = bestScoreD - qa.score();
                            qa.score(std::exp(-v));
                          });
          } else {
            numDecoyFrags += bestHitDecoy ? 1 : 0;
            ++numDropped;
            jointHitGroup.clearAlignments();
          }
        } //end-if validate mapping

        if (writeQuasimappings) {
          rapmap::utils::writeAlignmentsToStream(rp, formatter,
                                                 hctr, jointHits, sstream);
        }

      if (writeUnmapped and mapType == salmon::utils::MappingType::UNMAPPED) {
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

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<QuasiAlignment> hitLists = {structureVec.begin(), structureVec.begin()+rangeSize};
    processMiniBatchSimple<QuasiAlignment>(
                                     readExp, fmCalc, firstTimestepOfRound, rl, salmonOpts, hitLists,
                                     transcripts, clusterForest, fragLengthDist, observedBiasParams,
                                     numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac,
                                     alevinOpts, /*uniqueFLD,*/ iomutex);
  }

  if (maxZeroFrac > 5.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }

  mstats.numDecoyFragments += numDecoyFrags;
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
                        std::atomic<uint32_t>& smallSeqs,
                        std::atomic<uint32_t>& nSeqs,
                        bool initialRound,
                        std::atomic<bool>& burnedIn, ForgettingMassCalculator& fmCalc,
                        FragmentLengthDistribution& fragLengthDist,
                        SalmonOpts& salmonOpts, bool greedyChain,
                        std::mutex& iomutex, size_t numThreads,
                        std::vector<AlnGroupVec<AlnT>>& structureVec, volatile bool& writeToCache,
                        AlevinOpts<ProtocolT>& alevinOpts,
                        SoftMapT& barcodeMap,
                        spp::sparse_hash_map<std::string, uint32_t>& trBcs, MappingStatistics& mstats) {

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

    /*
    std::vector<std::vector<uint64_t>> uniqueFLDs(numThreads);
    for (size_t i = 0; i < numThreads; ++i) { uniqueFLDs[i] = std::vector<uint64_t>(100000, 0); }
    */

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
                                                                            upperBoundHits, smallSeqs, nSeqs, sidx->quasiIndexPerfectHash64(), transcripts,
                                                                            fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                                                                            salmonOpts, iomutex, initialRound,
                                                                            burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs,
                                                                            mstats/*, uniqueFLDs[i]*/);
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
                                                                          upperBoundHits, smallSeqs, nSeqs, sidx->quasiIndex64(), transcripts, fmCalc,
                                                                          clusterForest, fragLengthDist, observedBiasParams[i],
                                                                          salmonOpts, iomutex, initialRound,
                                                                          burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs,
                                                                          mstats/*, uniqueFLDs[i]*/);
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
                                                                            upperBoundHits, smallSeqs, nSeqs, sidx->quasiIndexPerfectHash32(), transcripts,
                                                                            fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                                                                            salmonOpts, iomutex, initialRound,
                                                                            burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs,
                                                                            mstats/*, uniqueFLDs[i]*/);
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
                                                                          upperBoundHits, smallSeqs, nSeqs, sidx->quasiIndex32(), transcripts, fmCalc,
                                                                          clusterForest, fragLengthDist, observedBiasParams[i],
                                                                          salmonOpts, iomutex, initialRound,
                                                                          burnedIn, writeToCache, alevinOpts, barcodeMap, trBcs,
                                                                          mstats/*, uniqueFLDs[i]*/);
          };
          threads.emplace_back(threadFun);
        }

      } // End spawn current thread

    } // End spawn all threads

    for (auto& t : threads) {
      t.join();
    }

    pairedParserPtr->stop();

    // At this point, if we were using decoy transcripts, we don't need them anymore and can get
    // rid of them.
    readExp.dropDecoyTranscripts();

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
                     MappingStatistics& mstats,
                     uint32_t numQuantThreads,
                     AlevinOpts<ProtocolT>& alevinOpts,
                     SoftMapT& barcodeMap,
                     spp::sparse_hash_map<std::string, uint32_t>& trBcs) {

  bool burnedIn = (salmonOpts.numBurninFrags == 0);
  uint64_t numRequiredFragments = salmonOpts.numRequiredFragments;
  std::atomic<uint64_t> upperBoundHits{0};
  auto& refs = experiment.transcripts();
  size_t numTranscripts = refs.size();
  // The *total* number of fragments observed so far (over all passes through
  // the data).
  std::atomic<uint64_t> numObservedFragments{0};
  std::atomic<uint32_t> smallSeqs{0};
  std::atomic<uint32_t> nSeqs{0};
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
                             upperBoundHits, smallSeqs, nSeqs, initialRound, burnedIn, fmCalc,
                             fragLengthDist, salmonOpts,
                             greedyChain, ioMutex,
                             numQuantThreads, groupVec, writeToCache,
                             alevinOpts, barcodeMap, trBcs, mstats);

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

  salmonOpts.jointLog->info("Number of fragments discarded because they are best-mapped to decoys : {:n}",
                            mstats.numDecoyFragments.load());

  if (totalAssignedFragments < salmonOpts.numBurninFrags) {
    std::atomic<bool> dummyBool{false};
    experiment.updateTranscriptLengthsAtomic(dummyBool);
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

  if (smallSeqs > 100) {
    jointLog->warn("Found {} reads with CB+UMI length smaller than expected.\n"
                   "Please report on github if this number is too large", smallSeqs);
  }
  if (nSeqs > 100) {
    jointLog->warn("Found {} reads with `N` in the UMI sequence and ignored the reads.\n"
                   "Please report on github if this number is too large", nSeqs);
  }
  alevinOpts.noisyUmis = nSeqs;
  alevinOpts.eqReads = totalAssignedFragments;
  alevinOpts.mappingRate = experiment.effectiveMappingRate() * 100.0;
  //+++++++++++++++++++++++++++++++++++++++
  jointLog->info("Mapping rate = {}\%\n",
                 experiment.effectiveMappingRate() * 100.0);
  jointLog->info("finished quantifyLibrary()");
}

template <typename ProtocolT>
void alevinOptimize( std::vector<std::string>& trueBarcodesVec,
                     spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                     spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                     EqMapT& fullEqMap,
                     AlevinOpts<ProtocolT>& aopt,
                     GZipWriter& gzw,
                     CFreqMapT& freqCounter,
                     size_t numLowConfidentBarcode) {
  std::vector<uint32_t> umiCount(trueBarcodesVec.size());
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

  ////////////////////////////////////////////
  // deduplication starts from here
  ////////////////////////////////////////////

  if(not aopt.noDedup) {
    aopt.jointLog->info("Starting optimizer\n\n");
    aopt.jointLog->flush();

    CollapsedCellOptimizer optimizer;
    bool optSuccess = optimizer.optimize(fullEqMap,
                                         txpToGeneMap,
                                         geneIdxMap,
                                         aopt,
                                         gzw,
                                         trueBarcodesVec,
                                         umiCount,
                                         freqCounter,
                                         numLowConfidentBarcode);
    if (!optSuccess) {
      aopt.jointLog->error(
                      "The optimization algorithm failed. This is likely the result of "
                      "bad input (or a bug). If you cannot track down the cause, please "
                      "report this issue on GitHub.");
      aopt.jointLog->flush();
      exit(74);
    }
    aopt.jointLog->info("Finished optimizer");
  }
  else{
    aopt.jointLog->warn("No Dedup command given, is it what you want?");
  }
}

template <typename ProtocolT>
int alevinQuant(AlevinOpts<ProtocolT>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
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

    MappingStatistics mstats;
    ReadExperimentT experiment(readLibraries, indexDirectory, sopt);
    //experiment.computePolyAPositions();

    // This will be the class in charge of maintaining our
    // rich equivalence classes
    experiment.equivalenceClassBuilder().setMaxResizeThreads(sopt.maxHashResizeThreads);
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

    std::vector<std::string> trueBarcodesVec (trueBarcodes.begin(),
                                              trueBarcodes.end());
    std::sort (trueBarcodesVec.begin(), trueBarcodesVec.end(),
               [&freqCounter, &jointLog](const std::string& i,
                                         const std::string& j){
                 uint32_t iCount, jCount;
                 auto itI = freqCounter.find(i);
                 auto itJ = freqCounter.find(j);
                 bool iOk = itI != freqCounter.end();
                 bool jOk = itJ != freqCounter.end();
                 if (not iOk or not jOk){
                   jointLog->error("Barcode not found in frequency table");
                   jointLog->flush();
                   exit(1);
                 }
                 iCount = *itI;
                 jCount = *itJ;
                 if (iCount > jCount){
                   return true;
                 }
                 else if (iCount < jCount){
                   return false;
                 }
                 else{
                   // stable sorting
                   if (i>j){
                     return true;
                   }
                   else{
                     return false;
                   }
                 }
               });
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
                                    mstats, sopt.numThreads, aopt,
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

    if(aopt.dumpBFH){
      gzw.writeBFH(aopt.outputDirectory, experiment,
                   aopt.protocol.umiLength, trueBarcodesVec);
    }

    if (aopt.dumpBarcodeEq and not aopt.noDedup){
      std::ofstream oFile;
      boost::filesystem::path oFilePath = aopt.outputDirectory / "cell_eq_order.txt";
      oFile.open(oFilePath.string());
      for (auto& bc : trueBarcodesVec) {
        oFile << bc << "\n";
      }
      oFile.close();

      {//dump transcripts names
        boost::filesystem::path tFilePath = aopt.outputDirectory / "transcripts.txt";
        std::ofstream tFile(tFilePath.string());
        for (auto& txp: experiment.transcripts()) {
          tFile << txp.RefName << "\n";
        }
        tFile.close();
      }
    }

    alevinOptimize(trueBarcodesVec, txpToGeneMap, geneIdxMap,
                   experiment.equivalenceClassBuilder().eqMap(),
                   aopt, gzw, freqCounter,
                   numLowConfidentBarcode);
    jointLog->flush();

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
    gzw.writeMeta(sopt, experiment, mstats);

    gzw.writeMetaAlevin(aopt, bfs::path(sopt.auxDir));

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
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::InDrop>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::ChromiumV3>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::Chromium>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::Gemcode>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::CELSeq>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::CELSeq2>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
template
int alevinQuant(AlevinOpts<apt::Custom>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode);
