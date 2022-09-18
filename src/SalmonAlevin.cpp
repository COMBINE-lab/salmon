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
#include <memory>
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
// #include "jellyfish/mer_dna.hpp"

// Boost Includes
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/filesystem.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/program_options.hpp>
#include <boost/range/irange.hpp>
#include <boost/thread/thread.hpp>

// TBB Includes
#include "oneapi/tbb/blocked_range.h"
#include "oneapi/tbb/concurrent_queue.h"
#include "oneapi/tbb/concurrent_unordered_map.h"
#include "oneapi/tbb/concurrent_unordered_set.h"
#include "oneapi/tbb/concurrent_vector.h"
#include "oneapi/tbb/parallel_for.h"
#include "oneapi/tbb/parallel_for_each.h"
#include "oneapi/tbb/parallel_reduce.h"
#include "oneapi/tbb/partitioner.h"

// logger includes
#include "spdlog/spdlog.h"

// Cereal includes
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"

#include "concurrentqueue.h"

#include <cuckoohash_map.hh>

// core includes
#include "core/range.hpp"

// radicl includes
#include "radicl/BasicBinWriter.hpp"
#include "radicl/RADHeader.hpp"

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
#include "SalmonMappingUtils.hpp"

#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"
#include "PairedAlignmentFormatter.hpp"

#include "pufferfish/Util.hpp"
#include "pufferfish/MemCollector.hpp"
#include "pufferfish/MemChainer.hpp"
#include "pufferfish/SAMWriter.hpp"
#include "pufferfish/PuffAligner.hpp"
#include "pufferfish/ksw2pp/KSW2Aligner.hpp"
#include "pufferfish/metro/metrohash64.h"
#include "parallel_hashmap/phmap.h"
#include "pufferfish/itlib/static_vector.hpp"
#include "pufferfish/SelectiveAlignmentUtils.hpp"

namespace alevin{

  /****** QUASI MAPPING DECLARATIONS *********/
  using MateStatus = pufferfish::util::MateStatus;
  using QuasiAlignment = pufferfish::util::QuasiAlignment;
  /****** QUASI MAPPING DECLARATIONS  *******/

  using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
  using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

  using TranscriptID = uint32_t;
  using TranscriptIDVector = std::vector<TranscriptID>;
  using KmerIDMap = std::vector<TranscriptIDVector>;

  constexpr uint32_t miniBatchSize{5000};

  // NOTE: Consider the implications of using uint32_t below.
  // Currently, this should not limit the barcode length to <= 16 in 
  // "fry" mode (either sla or sketch) since we never actually put the 
  // barcode in the corresponding variable of the AlignmentGroup class.
  // In traditional alevin, this should be referring to a barcode *index*
  // rather than the sequence itself, so it should be OK so long as we have
  // < std::numeric_limits<uint32_t>::max() distinct barcodes / reads. 
  // However, that assumption should be tested more thoroughly.
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
  using AlnGroupQueue = oneapi::tbb::concurrent_queue<AlevinAlnGroup<AlnT>*>;
#endif

  //#include "LightweightAlignmentDefs.hpp"
}

//have to create new namespace because of multiple definition
using namespace alevin;

/* ALEVIN DECLERATIONS*/
using bcEnd = BarcodeEnd;
namespace aut = alevin::utils;
using BlockedIndexRange = oneapi::tbb::blocked_range<size_t>;
using ReadExperimentT = ReadExperiment<EquivalenceClassBuilder<SCTGValue>>;
/////// REDUNDANT CODE END//

template <typename AlnT>
void processMiniBatchSimple(ReadExperimentT& readExp, ForgettingMassCalculator& fmCalc,
                      ReadLibrary& readLib,
                      const SalmonOpts& salmonOpts,
                      AlnGroupVecRange<AlnT> batchHits,
                      std::vector<Transcript>& transcripts,
                      ClusterForest& clusterForest,
                      FragmentLengthDistribution& fragLengthDist,
                      std::atomic<uint64_t>& numAssignedFragments,
                      bool initialRound,
                      std::atomic<bool>& burnedIn, double& maxZeroFrac){

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

/**
 * Perform sketch mapping of sc reads and produce output file.
 **/
template <typename IndexT, typename ProtocolT>
void process_reads_sc_sketch(paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
                       AlnGroupVec<QuasiAlignment>& structureVec,
                       std::atomic<uint64_t>& numObservedFragments,
                       std::atomic<uint64_t>& numAssignedFragments,
                       std::atomic<uint64_t>& numUniqueMappings,
                       std::atomic<uint64_t>& validHits, 
                       std::atomic<uint32_t>& smallSeqs,
                       std::atomic<uint32_t>& nSeqs,
                       IndexT* qidx, std::vector<Transcript>& transcripts,
                       FragmentLengthDistribution& fragLengthDist, 
                       SalmonOpts& salmonOpts,
                       std::atomic<uint64_t>& num_chunks,
                       std::ofstream& rad_file,
                       std::ofstream& ubc_file, // unmapped barcode count
                       std::mutex& fileMutex,
                       std::mutex& ubc_file_mutex, // mutex for barcode count
                       std::mutex& iomutex, 
                       AlevinOpts<ProtocolT>& alevinOpts,
                       MappingStatistics& mstats) {
  (void) numAssignedFragments;
  (void) fragLengthDist;
  uint64_t count_fwd = 0, count_bwd = 0;
  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};
  salmon::utils::ShortFragStats shortFragStats;
  double maxZeroFrac{0.0};

  bool quiet = salmonOpts.quiet;
  if (salmonOpts.writeQualities) {
    salmonOpts.jointLog->warn("The --writeQualities flag has no effect when writing results to a RAD file.");
  }

  BasicBinWriter bw;
  uint32_t num_reads_in_chunk{0};
  bw << num_reads_in_chunk;
  bw << num_reads_in_chunk;

  size_t minK = qidx->k();
  int32_t signed_k = static_cast<int32_t>(minK);
  size_t locRead{0};
  size_t rangeSize{0};
  uint64_t localNumAssignedFragments{0};
  uint64_t localNumUniqueMappings{0};

  size_t readLenLeft{0};
  size_t readLenRight{0};

  salmon::utils::MappingType mapType{salmon::utils::MappingType::UNMAPPED};

  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);
  PairedAlignmentFormatter<IndexT*> formatter(qidx);
  
  // map to recall the number of unmapped reads we see 
  // for each barcode
  phmap::flat_hash_map<uint64_t, uint32_t> unmapped_bc_map; 

  // check the frequency and decide here if we should
  // be attempting recovery of highly-multimapping reads
  const size_t max_occ_default = salmonOpts.maxReadOccs;
  const size_t max_occ_recover = salmonOpts.maxRecoverReadOccs;
  const bool attempt_occ_recover = (max_occ_recover > max_occ_default);

  // think if this should be different since it's at the read level
  const size_t max_read_occs = salmonOpts.maxReadOccs;
  size_t alt_max_occ = max_read_occs;

  // we'll use this for mapping
  salmon::mapping_utils::pasc::mapping_cache_info<IndexT> map_cache_first(qidx, max_occ_default, max_occ_recover, max_read_occs);
  salmon::mapping_utils::pasc::mapping_cache_info<IndexT> map_cache_second(qidx, max_occ_default, max_occ_recover, max_read_occs);

  std::vector<salmon::mapping_utils::pasc::simple_hit>* accepted_hits_left_ptr;
  std::vector<salmon::mapping_utils::pasc::simple_hit>* accepted_hits_right_ptr;
  std::vector<salmon::mapping_utils::pasc::simple_hit>* accepted_hits_ptr;

  // buffer in which to accumulate accepted hits if we are merging for left 
  // and right reads
  std::vector<salmon::mapping_utils::pasc::simple_hit> accepted_hits_buffer;

  std::string readBuffer;
  std::string umi(alevinOpts.protocol.umiLength, 'N');
  std::string barcode(alevinOpts.protocol.barcodeLength, 'N');
  //////////////////////

  auto localProtocol = alevinOpts.protocol;
  // start parsing reads
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

    LibraryFormat expectedLibraryFormat = rl.format();

    std::string extraBAMtags;
    if(writeQuasimappings) {
    	size_t reserveSize { alevinOpts.protocol.barcodeLength + alevinOpts.protocol.umiLength + 12};
    	extraBAMtags.reserve(reserveSize);
    }

    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      auto& rp = rg[i];
      readLenLeft = rp.first.seq.length();
      readLenRight= rp.second.seq.length();

      // while the pasc code path will use the jointHitGroup to 
      // store UMI and BarCode info, it won't actually use the 
      // "alignment" fields
      bool tooShortRight = (readLenRight < (minK+alevinOpts.trimRight));
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();

      mapType = salmon::utils::MappingType::UNMAPPED;
      // ensure that accepted_hits_ptr points to a valid 
      // vector in case it's not set below
      accepted_hits_ptr = &accepted_hits_buffer;
      alt_max_occ = max_read_occs;

      //////////////////////////////////////////////////////////////
      // extracting barcodes
      size_t barcodeLength = alevinOpts.protocol.barcodeLength;
      size_t umiLength = alevinOpts.protocol.umiLength;

      nonstd::optional<uint32_t> barcodeIdx;
      extraBAMtags.clear();
      bool seqOk = false;

      if (alevinOpts.protocol.end == bcEnd::FIVE ||
          alevinOpts.protocol.end == bcEnd::THREE) {
        // If the barcode sequence could be extracted, then this is set to true,
        // but the barcode sequence itself may still be invalid (e.g. contain `N` characters).
        // However, if extracted_barcode is false here, there is no hope to even recover the
        // barcode and we shouldn't attempt it.
        bool extracted_barcode = aut::extractBarcode(rp.first.seq, rp.second.seq, localProtocol, barcode);
        // If we could pull out something where the barcode sequence should have been
        // then continue to process it.
        if (extracted_barcode) {
          // if the barcode consisted of valid nucleotides, then seqOk is true
          // otherwise false
          seqOk =  aut::sequenceCheck(barcode, Sequence::BARCODE);
          if (not seqOk){
            // If the barcode contained invalid nucleotides
            // this attempts to replace the first one with an `A`.
            // If this returns true, there was only one `N` and we
            // replaced it; otherwise there was more than one `N`
            // and the barcode sequence should be treated as invalid.
            seqOk = aut::recoverBarcode(barcode);
          }
        }

        // If we have a valid barcode
        if (seqOk) {
          bool umi_ok = aut::extractUMI(rp.first.seq, rp.second.seq, localProtocol, umi);
          if ( !umi_ok ) {
            smallSeqs += 1;
          } else {
            alevin::types::AlevinUMIKmer umiIdx;
            bool isUmiIdxOk = umiIdx.fromChars(umi);

            if (isUmiIdxOk) {
              jointHitGroup.setUMI(umiIdx.word(0));

              // get the subsequence for the alignable portion of the read
              std::string* readSubSeq = aut::getReadSequence(localProtocol, rp.first.seq, rp.second.seq, readBuffer);

              // try and map it
              bool early_exit = salmon::mapping_utils::pasc::map_read(readSubSeq, map_cache_second, PairingStatus::UNPAIRED_RIGHT);
              // in this case accepted_hits_ptr gets set to what we found
              accepted_hits_ptr = &map_cache_second.accepted_hits;
              alt_max_occ = map_cache_second.alt_max_occ;
            } else {
              nSeqs += 1;
            }
          }
        }
      } else {
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

      // If the read mapped to > maxReadOccs places, discard it
      if (accepted_hits_ptr->size() > alt_max_occ) { 
        accepted_hits_ptr->clear();
      } else if (!accepted_hits_ptr->empty()) {
        mapType = salmon::utils::MappingType::SINGLE_MAPPED;
      }


      alevin::types::AlevinCellBarcodeKmer bck;
      bool barcode_ok = bck.fromChars(barcode);
 
      // NOTE: Think if we should put decoy mappings in the RAD file
      if (mapType == salmon::utils::MappingType::SINGLE_MAPPED) {
        // num aln
        bw << static_cast<uint32_t>(accepted_hits_ptr->size()); 

        // bc
        // if we can fit the barcode into an integer 
        if ( barcodeLength <= 32 ) { 
          if (barcode_ok) {
            if ( barcodeLength <= 16 ) { // can use 32-bit int
              uint32_t shortbck = static_cast<uint32_t>(0x00000000FFFFFFFF & bck.word(0));
              bw << shortbck;
            } else { // must use 64-bit int
              bw << bck.word(0);
            }
          }
        } else { // must use a string for the barcode
          bw << barcode;
        }

        // umi
        if ( umiLength <= 16 ) { // if we can use 32-bit int 
          uint64_t umiint = jointHitGroup.umi();
          uint32_t shortumi = static_cast<uint32_t>(0x00000000FFFFFFFF & umiint);
          bw << shortumi;
        } else if ( umiLength <= 32 ) { // if we can use 64-bit int
          uint64_t umiint = jointHitGroup.umi();
          bw << umiint;
        } else { // must use string
          bw << umi;
        }

        for (auto& aln : (*accepted_hits_ptr)) {
          uint32_t fw_mask = aln.is_fw ? 0x80000000 : 0x00000000;
          bw << (aln.tid | fw_mask);
        }
        ++num_reads_in_chunk;
      } else { // if read was not mapped
        if (barcode_ok) {
          unmapped_bc_map[bck.word(0)] += 1;
        }
      }


      if (num_reads_in_chunk > 5000) {
        ++num_chunks;
        uint32_t num_bytes = bw.num_bytes();
        bw.write_integer_at_offset(0, num_bytes);
        bw.write_integer_at_offset(sizeof(num_bytes), num_reads_in_chunk);
        fileMutex.lock();
        rad_file << bw;
        fileMutex.unlock();
        bw.clear();
        num_reads_in_chunk = 0;

        // reserve space for headers of next chunk
        bw << num_reads_in_chunk;
        bw << num_reads_in_chunk;
      }

      validHits += accepted_hits_ptr->size();
      localNumAssignedFragments += (accepted_hits_ptr->size() > 0);
      localNumUniqueMappings += (accepted_hits_ptr->size() == 1) ? 1 : 0;
      locRead++;
      ++numObservedFragments;
      jointHitGroup.clearAlignments();

      if (!quiet and numObservedFragments % 500000 == 0) {
        iomutex.lock();
        const char RESET_COLOR[] = "\x1b[0m";
        char green[] = "\x1b[30m";
        green[3] = '0' + static_cast<char>(fmt::GREEN);
        char red[] = "\x1b[30m";
        red[3] = '0' + static_cast<char>(fmt::RED);
        fmt::print(
            stderr, "\033[A\r\r{}processed{} {} Million {}fragments{}\n",
            green, red, numObservedFragments / 1000000, green, RESET_COLOR);
        fmt::print(stderr, "hits: {}, hits per frag:  {}", validHits,
                   validHits / static_cast<float>(prevObservedFrags));
        iomutex.unlock();
      }

    } // end for i < j->nb_filled

    prevObservedFrags = numObservedFragments;
  }

  // If we have any unwritten records left to write
  if (num_reads_in_chunk > 0) {
    ++num_chunks;
    uint32_t num_bytes = bw.num_bytes();
    bw.write_integer_at_offset(0, num_bytes);
    bw.write_integer_at_offset(sizeof(num_bytes), num_reads_in_chunk);
    fileMutex.lock();
    rad_file << bw;
    fileMutex.unlock();
    bw.clear();
    num_reads_in_chunk = 0;
  }

  // unmapped barcode writer
  { // make a scope and dump the unmapped barcode counts
    BasicBinWriter ubcw;
    for (auto& kv : unmapped_bc_map) {
      ubcw << kv.first;
      ubcw << kv.second;
    }
    ubc_file_mutex.lock();
    ubc_file << ubcw;
    ubc_file_mutex.unlock();
    ubcw.clear();
  }

  if (maxZeroFrac > 5.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }
  numAssignedFragments += localNumAssignedFragments;
  numUniqueMappings += localNumUniqueMappings;
  mstats.numDecoyFragments += 0; // no decoys with pasc (yet?)
  readExp.updateShortFrags(shortFragStats);
}

/**
 * Perform selective alignment of sc reads and produce output file.
 **/
template <typename IndexT, typename ProtocolT>
void process_reads_sc_align(paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
                       AlnGroupVec<QuasiAlignment>& structureVec,
                       std::atomic<uint64_t>& numObservedFragments,
                       std::atomic<uint64_t>& numAssignedFragments,
                       std::atomic<uint64_t>& validHits, 
                       std::atomic<uint32_t>& smallSeqs,
                       std::atomic<uint32_t>& nSeqs,
                       IndexT* qidx, std::vector<Transcript>& transcripts,
                       FragmentLengthDistribution& fragLengthDist, 
                       SalmonOpts& salmonOpts,
                       std::atomic<uint64_t>& num_chunks,
                       std::ofstream& rad_file,
                       std::ofstream& ubc_file, // unmapped barcode count
                       std::mutex& fileMutex,
                       std::mutex& ubc_file_mutex, // mutex for barcode count
                       std::mutex& iomutex, 
                       AlevinOpts<ProtocolT>& alevinOpts,
                       MappingStatistics& mstats) {
  (void) numAssignedFragments;
  (void) fragLengthDist;
  uint64_t count_fwd = 0, count_bwd = 0;
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

  BasicBinWriter bw;
  uint32_t num_reads_in_chunk{0};
  bw << num_reads_in_chunk;
  bw << num_reads_in_chunk;

  size_t minK = qidx->k();
  size_t locRead{0};
  size_t rangeSize{0};
  uint64_t localNumAssignedFragments{0};

  bool consistentHits = salmonOpts.consistentHits;
  bool quiet = salmonOpts.quiet;

  size_t maxNumHits{salmonOpts.maxReadOccs};
  size_t readLenLeft{0};
  size_t readLenRight{0};

  constexpr const int32_t invalidScore = std::numeric_limits<int32_t>::min();
  MemCollector<IndexT> memCollector(qidx);
  ksw2pp::KSW2Aligner aligner;
  pufferfish::util::AlignmentConfig aconf;
  pufferfish::util::MappingConstraintPolicy mpol;
  bool initOK = salmon::mapping_utils::initMapperSettings(salmonOpts, memCollector, aligner, aconf, mpol);
  PuffAligner puffaligner(qidx->refseq_, qidx->refAccumLengths_, qidx->k(), aconf, aligner);

  pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> hits;
  std::vector<pufferfish::util::MemCluster> recoveredHits;
  std::vector<pufferfish::util::JointMems> jointHits;
  PairedAlignmentFormatter<IndexT*> formatter(qidx);
  if (salmonOpts.writeQualities) {
    salmonOpts.jointLog->warn("The --writeQualities flag has no effect when writing results to a RAD file.");
  }
  pufferfish::util::QueryCache qc;

  bool mimicStrictBT2 = salmonOpts.mimicStrictBT2;
  bool mimicBT2 = salmonOpts.mimicBT2;
  bool noDovetail = !salmonOpts.allowDovetail;
  bool useChainingHeuristic = !salmonOpts.disableChainingHeuristic;

  pufferfish::util::HitCounters hctr;
  salmon::utils::MappingType mapType{salmon::utils::MappingType::UNMAPPED};
  bool hardFilter = salmonOpts.hardFilter;

  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);

  // map to recall the number of unmapped reads we see 
  // for each barcode
  phmap::flat_hash_map<uint64_t, uint32_t> unmapped_bc_map; 

  //////////////////////
  // NOTE: validation mapping based new parameters
  std::string rc1; rc1.reserve(300);
  AlnCacheMap alnCache; alnCache.reserve(16);
  /*
  auto ap{selective_alignment::utils::AlignmentPolicy::DEFAULT};
  if (mimicBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2;
  } else if (mimicStrictBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2_STRICT;
  }
  */

  size_t numMappingsDropped{0};
  size_t numDecoyFrags{0};
  const double decoyThreshold = salmonOpts.decoyThreshold;

  salmon::mapping_utils::MappingScoreInfo msi(decoyThreshold);
  // we only collect detailed decoy information if we will be 
  // writing output to SAM.
  msi.collect_decoys(writeQuasimappings);

  std::string* readSubSeq{nullptr};
  std::string readBuffer;
  std::string umi(alevinOpts.protocol.umiLength, 'N');
  std::string barcode(alevinOpts.protocol.barcodeLength, 'N');
  //////////////////////

  bool tryAlign{salmonOpts.validateMappings};
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

    LibraryFormat expectedLibraryFormat = rl.format();

    std::string extraBAMtags;
    if(writeQuasimappings) {
    	size_t reserveSize { alevinOpts.protocol.barcodeLength + alevinOpts.protocol.umiLength + 12};
    	extraBAMtags.reserve(reserveSize);
    }

    auto localProtocol = alevinOpts.protocol;
    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      auto& rp = rg[i];
      readLenLeft = rp.first.seq.length();
      readLenRight= rp.second.seq.length();

      bool tooShortRight = (readLenRight < (minK+alevinOpts.trimRight));
      //localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();
      auto& jointAlignments= jointHitGroup.alignments();

      hits.clear();
      jointHits.clear();
      memCollector.clear();
      //jointAlignments.clear();
      readSubSeq = nullptr;//.clear();
      mapType = salmon::utils::MappingType::UNMAPPED;

      //////////////////////////////////////////////////////////////
      // extracting barcodes
      size_t barcodeLength = alevinOpts.protocol.barcodeLength;
      size_t umiLength = alevinOpts.protocol.umiLength;
      //umi.clear();
      //barcode.clear();
      nonstd::optional<uint32_t> barcodeIdx;
      extraBAMtags.clear();
      bool seqOk = false;

      if (alevinOpts.protocol.end == bcEnd::FIVE ||
          alevinOpts.protocol.end == bcEnd::THREE){
        // If the barcode sequence could be extracted, then this is set to true,
        // but the barcode sequence itself may still be invalid (e.g. contain `N` characters).
        // However, if extracted_barcode is false here, there is no hope to even recover the
        // barcode and we shouldn't attempt it.
        bool extracted_barcode = aut::extractBarcode(rp.first.seq, rp.second.seq, localProtocol, barcode);
        // If we could pull out something where the barcode sequence should have been
        // then continue to process it.
        if (extracted_barcode) {
          // if the barcode consisted of valid nucleotides, then seqOk is true
          // otherwise false
          seqOk =  aut::sequenceCheck(barcode, Sequence::BARCODE);
          if (not seqOk){
            // If the barcode contained invalid nucleotides
            // this attempts to replace the first one with an `A`.
            // If this returns true, there was only one `N` and we
            // replaced it; otherwise there was more than one `N`
            // and the barcode sequence should be treated as invalid.
            seqOk = aut::recoverBarcode(barcode);
          }
        }

        // If we have a valid barcode
        if (seqOk) {
          bool umi_ok = aut::extractUMI(rp.first.seq, rp.second.seq, localProtocol, umi);
          if ( !umi_ok ) {
            smallSeqs += 1;
          } else {
            alevin::types::AlevinUMIKmer umiIdx;
            bool isUmiIdxOk = umiIdx.fromChars(umi);

            if (isUmiIdxOk) {
              jointHitGroup.setUMI(umiIdx.word(0));
              /*
              auto seq_len = rp.second.seq.size();
              if (alevinOpts.trimRight > 0) {
                if (!tooShortRight) {
                  readSubSeq = rp.second.seq.substr(0, seq_len - alevinOpts.trimRight);
                  auto rh = memCollector(readSubSeq, qc,
                                         true, // isLeft
                                         false // verbose
                  );
                }
              } else {
                */
              readSubSeq = aut::getReadSequence(localProtocol, rp.first.seq, rp.second.seq, readBuffer);
              auto rh = tooShortRight ? false
                                      : memCollector(*readSubSeq, qc,
                                                     true, // isLeft
                                                     false // verbose
                                        );
              // }
              memCollector.findChains(
                  *readSubSeq, hits, salmonOpts.fragLenDistMax,
                  MateStatus::PAIRED_END_RIGHT,
                  useChainingHeuristic, // heuristic chaining
                  true,                 // isLeft
                  false                 // verbose
              );

              pufferfish::util::joinReadsAndFilterSingle(
                  hits, jointHits, readSubSeq->length(),
                  memCollector.getConsensusFraction());
            } else {
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


      // If the read mapped to > maxReadOccs places, discard it
      if (jointHits.size() > salmonOpts.maxReadOccs) {
        jointHitGroup.clearAlignments();
      }

        if (!jointHits.empty()) {
          puffaligner.clear();
          puffaligner.getScoreStatus().reset();
          msi.clear(jointHits.size());

          size_t idx{0};
          bool is_multimapping = (jointHits.size() > 1);

          for (auto &&jointHit : jointHits) {
            // for alevin, currently, we need these to have a mate status of PAIRED_END_RIGHT
            jointHit.mateStatus = MateStatus::PAIRED_END_RIGHT;
            auto hitScore = puffaligner.calculateAlignments(*readSubSeq, jointHit, hctr, is_multimapping, salmonOpts.jointLog, false);
            bool validScore = (hitScore != invalidScore);
            numMappingsDropped += validScore ? 0 : 1;
            auto tid = qidx->getRefId(jointHit.tid);
            
            // NOTE: Here, we know that the read arising from the transcriptome is the "right"
            // read (read 2).  So we interpret compatibility in that context.
            // TODO: Make this code more generic and modular (account for the possibility of different library)
            // protocols or setups where the reads are not always "paired-end" and the transcriptomic read is not
            // always read 2 (@k3yavi).
            bool isCompat = (expectedLibraryFormat.strandedness == ReadStrandedness::U) or 
                            (jointHit.orphanClust()->isFw and (expectedLibraryFormat.strandedness == ReadStrandedness::AS)) or
                            (!jointHit.orphanClust()->isFw and (expectedLibraryFormat.strandedness == ReadStrandedness::SA));

            salmon::mapping_utils::updateRefMappings(tid, hitScore, isCompat, idx, transcripts, invalidScore, 
                                                     msi);
            ++idx;
          }

          bool bestHitDecoy = msi.haveOnlyDecoyMappings();
          if (msi.bestScore > invalidScore and !bestHitDecoy) {
            salmon::mapping_utils::filterAndCollectAlignments(jointHits,
                                                              readSubSeq->length(),
                                                              readSubSeq->length(),
                                                              false, // true for single-end false otherwise
                                                              tryAlign,
                                                              hardFilter,
                                                              salmonOpts.scoreExp,
                                                              salmonOpts.minAlnProb,
                                                              msi,
                                                              jointAlignments);
            if (!jointAlignments.empty()) {
              mapType = salmon::utils::MappingType::SINGLE_MAPPED;
            }
          } else {
            numDecoyFrags += bestHitDecoy ? 1 : 0;
            mapType = (bestHitDecoy) ? salmon::utils::MappingType::DECOY : salmon::utils::MappingType::UNMAPPED;
            if (bestHitDecoy) {
              salmon::mapping_utils::filterAndCollectAlignmentsDecoy(
                  jointHits, readSubSeq->length(),
                  readSubSeq->length(),
                  false, // true for single-end false otherwise
                  tryAlign, hardFilter, salmonOpts.scoreExp,
                  salmonOpts.minAlnProb, msi,
                  jointAlignments);
            } else {
              jointHitGroup.clearAlignments();
            }
          }
        } //end-if validate mapping

        alevin::types::AlevinCellBarcodeKmer bck;
        bool barcode_ok = bck.fromChars(barcode);

        // NOTE: Think if we should put decoy mappings in the RAD file
        if (mapType == salmon::utils::MappingType::SINGLE_MAPPED) {
          // na
          bw << static_cast<uint32_t>(jointAlignments.size()); 
          // read length
          // bw << static_cast<uint16_t>(readLenRight); 

          // bc
          // if we can if the barcode into an integer 
          if ( barcodeLength <= 32 ) { 

            if (barcode_ok) {
              //alevinOpts.jointLog->info("BARCODE : {} \t ENC : {} ", barcode, bck.word(0));
              if ( barcodeLength <= 16 ) { // can use 32-bit int
                uint32_t shortbck = static_cast<uint32_t>(0x00000000FFFFFFFF & bck.word(0));
                //alevinOpts.jointLog->info("shortbck : {} ", shortbck);
                bw << shortbck;
              } else { // must use 64-bit int
                bw << bck.word(0);
              }
            }
          } else { // must use a string for the barcode
            bw << barcode;
          }

          // umi
          if ( umiLength <= 16 ) { // if we can use 32-bit int 
            uint64_t umiint = jointHitGroup.umi();
            uint32_t shortumi = static_cast<uint32_t>(0x00000000FFFFFFFF & umiint);
            bw << shortumi;
          } else if ( umiLength <= 32 ) { // if we can use 64-bit int
            uint64_t umiint = jointHitGroup.umi();
            bw << umiint;
          } else { // must use string
            bw << umi;
          }


          for (auto& aln : jointAlignments) {
            uint32_t fw_mask = aln.fwd ? 0x80000000 : 0x00000000;
            //bw << is_fw;
            bw << (aln.tid | fw_mask);
            // position 
            //bw << static_cast<uint32_t>((aln.pos < 0) ? 0 : aln.pos);
          }
          ++num_reads_in_chunk;
        } else {
          if (barcode_ok) {
            unmapped_bc_map[bck.word(0)] += 1;
          } 
        }


        if (num_reads_in_chunk > 5000) {
          ++num_chunks;
          uint32_t num_bytes = bw.num_bytes();
          bw.write_integer_at_offset(0, num_bytes);
          bw.write_integer_at_offset(sizeof(num_bytes), num_reads_in_chunk);
          fileMutex.lock();
          rad_file << bw;
          fileMutex.unlock();
          bw.clear();
          num_reads_in_chunk = 0;

          // reserve space for headers of next chunk
          bw << num_reads_in_chunk;
          bw << num_reads_in_chunk;
        }

        validHits += jointAlignments.size();
        localNumAssignedFragments += (jointAlignments.size() > 0);
        locRead++;
        ++numObservedFragments;
        jointHitGroup.clearAlignments();

        if (!quiet and numObservedFragments % 500000 == 0) {
          iomutex.lock();
          const char RESET_COLOR[] = "\x1b[0m";
          char green[] = "\x1b[30m";
          green[3] = '0' + static_cast<char>(fmt::GREEN);
          char red[] = "\x1b[30m";
          red[3] = '0' + static_cast<char>(fmt::RED);
          fmt::print(
              stderr, "\033[A\r\r{}processed{} {} Million {}fragments{}\n",
              green, red, numObservedFragments / 1000000, green, RESET_COLOR);
          fmt::print(stderr, "hits: {}, hits per frag:  {}", validHits,
                     validHits / static_cast<float>(prevObservedFrags));
          iomutex.unlock();
        }

    } // end for i < j->nb_filled

    prevObservedFrags = numObservedFragments;
    /*
    AlnGroupVecRange<QuasiAlignment> hitLists = {structureVec.begin(), structureVec.begin()+rangeSize};
    processMiniBatchSimple<QuasiAlignment>(
                                     readExp, fmCalc, rl, salmonOpts, hitLists,
                                     transcripts, clusterForest, fragLengthDist,
                                     numAssignedFragments, initialRound, burnedIn, maxZeroFrac);
                                     */
  }

  // If we have any unwritten records left to write
  if (num_reads_in_chunk > 0) {
    ++num_chunks;
    uint32_t num_bytes = bw.num_bytes();
    bw.write_integer_at_offset(0, num_bytes);
    bw.write_integer_at_offset(sizeof(num_bytes), num_reads_in_chunk);
    fileMutex.lock();
    rad_file << bw;
    fileMutex.unlock();
    bw.clear();
    num_reads_in_chunk = 0;
  }

    
  // unmapped barcode writer
  { // make a scope and dump the unmapped barcode counts
    BasicBinWriter ubcw;
    for (auto& kv : unmapped_bc_map) {
      ubcw << kv.first;
      ubcw << kv.second;
    }
    ubc_file_mutex.lock();
    ubc_file << ubcw;
    ubc_file_mutex.unlock();
    ubcw.clear();
  }

  if (maxZeroFrac > 5.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }
  numAssignedFragments += localNumAssignedFragments;
  mstats.numDecoyFragments += numDecoyFrags;
  readExp.updateShortFrags(shortFragStats);
}
// END SC_ALIGN

/// START QUASI
template <typename IndexT, typename ProtocolT>
void processReadsQuasi(
                       paired_parser* parser, ReadExperimentT& readExp, ReadLibrary& rl,
                       AlnGroupVec<QuasiAlignment>& structureVec,
                       std::atomic<uint64_t>& numObservedFragments,
                       std::atomic<uint64_t>& numAssignedFragments,
                       std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
                       std::atomic<uint32_t>& smallSeqs,
                       std::atomic<uint32_t>& nSeqs,
                       IndexT* qidx, std::vector<Transcript>& transcripts,
                       ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
                       FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
                       SalmonOpts& salmonOpts,
                       std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
                       AlevinOpts<ProtocolT>& alevinOpts,
                       SoftMapT& barcodeMap,
                       spp::sparse_hash_map<std::string, uint32_t>& trBcs,
                       MappingStatistics& mstats
                       /*,std::vector<uint64_t>& uniqueFLD*/) {
  uint64_t count_fwd = 0, count_bwd = 0;
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
  // Mer leftMer;
  // Mer rightMer;

  //auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
  size_t minK = qidx->k();

  size_t locRead{0};
  //uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};
  uint64_t localNumAssignedFragments{0};
  bool consistentHits = salmonOpts.consistentHits;
  bool quiet = salmonOpts.quiet;

  size_t maxNumHits{salmonOpts.maxReadOccs};
  size_t readLenLeft{0};
  size_t readLenRight{0};

  constexpr const int32_t invalidScore = std::numeric_limits<int32_t>::min();
  MemCollector<IndexT> memCollector(qidx);
  ksw2pp::KSW2Aligner aligner;
  pufferfish::util::AlignmentConfig aconf;
  pufferfish::util::MappingConstraintPolicy mpol;
  bool initOK = salmon::mapping_utils::initMapperSettings(salmonOpts, memCollector, aligner, aconf, mpol);
  PuffAligner puffaligner(qidx->refseq_, qidx->refAccumLengths_, qidx->k(), aconf, aligner);

  pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> hits;
  std::vector<pufferfish::util::MemCluster> recoveredHits;
  std::vector<pufferfish::util::JointMems> jointHits;
  PairedAlignmentFormatter<IndexT*> formatter(qidx);

  // Says if we should check that quality values exist
  // in the case the user requested to `--writeQualities`,
  // because they may have accidentially passed in a FASTA
  // file.
  bool check_qualities = true;
  if (salmonOpts.writeQualities) {
    formatter.enable_qualities();
  } else {
    // we don't have to worry about
    // checking qualities because 
    // we aren't writing them.
    check_qualities = false;
    formatter.disable_qualities();
  } 
  pufferfish::util::QueryCache qc;

  bool mimicStrictBT2 = salmonOpts.mimicStrictBT2;
  bool mimicBT2 = salmonOpts.mimicBT2;
  bool noDovetail = !salmonOpts.allowDovetail;
  bool useChainingHeuristic = !salmonOpts.disableChainingHeuristic;

  pufferfish::util::HitCounters hctr;
  salmon::utils::MappingType mapType{salmon::utils::MappingType::UNMAPPED};
  bool hardFilter = salmonOpts.hardFilter;

  fmt::MemoryWriter sstream;
  auto* qmLog = salmonOpts.qmLog.get();
  bool writeQuasimappings = (qmLog != nullptr);
  // if we aren't writing output at all, don't bother
  // checking for quality scores either.
  if (!writeQuasimappings) { check_qualities = false; }

  //////////////////////
  // NOTE: validation mapping based new parameters
  std::string rc1; rc1.reserve(300);
  AlnCacheMap alnCache; alnCache.reserve(16);

  /*
  auto ap{selective_alignment::utils::AlignmentPolicy::DEFAULT};
  if (mimicBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2;
  } else if (mimicStrictBT2) {
    ap = selective_alignment::utils::AlignmentPolicy::BT2_STRICT;
  }
  */

  size_t numMappingsDropped{0};
  size_t numDecoyFrags{0};
  const double decoyThreshold = salmonOpts.decoyThreshold;

  salmon::mapping_utils::MappingScoreInfo msi(decoyThreshold);
  // we only collect detailed decoy information if we will be 
  // writing output to SAM.
  msi.collect_decoys(writeQuasimappings);

  std::string* readSubSeq{nullptr};
  std::string readBuffer;
  std::string umi(alevinOpts.protocol.umiLength, 'N');
  std::string barcode(alevinOpts.protocol.barcodeLength, 'N');
  //////////////////////

  bool tryAlign{salmonOpts.validateMappings};
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

    LibraryFormat expectedLibraryFormat = rl.format();

    std::string extraBAMtags;
    if(writeQuasimappings) {
    	size_t reserveSize { alevinOpts.protocol.barcodeLength + alevinOpts.protocol.umiLength + 12};
    	extraBAMtags.reserve(reserveSize);
    }

    // if we need to disable writing quality values 
    // because the user passed in a FASTA file, do that
    // check here.
    if (check_qualities and (rangeSize > 0)) {
      auto& rp = rg[0];
      // a valid FASTQ record can't have an 
      // empty quality string, so then we will
      // treat this as a FASTA.
      if (rp.first.qual.empty() or rp.second.qual.empty()) { 
        formatter.disable_qualities();
        salmonOpts.jointLog->warn("The flag --writeQualities was provided,\n"
        "but read records (e.g. {}/{}) appear not to have quality strings!\n"
        "The input is being interpreted as a FASTA file, and the writing\n"
        "of quality scores is being disabled.\n", rp.first.name, rp.second.name);
      }
      // we won't bother to perform this check more than once.
      check_qualities = false;
    }

    auto localProtocol = alevinOpts.protocol;
    for (size_t i = 0; i < rangeSize; ++i) { // For all the read in this batch
      auto& rp = rg[i];
      readLenLeft = rp.first.seq.length();
      readLenRight= rp.second.seq.length();

      bool tooShortRight = (readLenRight < (minK+alevinOpts.trimRight));
      //localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();
      auto& jointAlignments= jointHitGroup.alignments();

      hits.clear();
      jointHits.clear();
      memCollector.clear();
      //jointAlignments.clear();
      //readSubSeq.clear();
      readSubSeq = nullptr;
      mapType = salmon::utils::MappingType::UNMAPPED;

      //////////////////////////////////////////////////////////////
      // extracting barcodes
      size_t barcodeLength = alevinOpts.protocol.barcodeLength;
      size_t umiLength = alevinOpts.protocol.umiLength;
      //umi.clear();
      //barcode.clear();
      nonstd::optional<uint32_t> barcodeIdx;
      extraBAMtags.clear();
      bool seqOk = false;

      if (alevinOpts.protocol.end == bcEnd::FIVE ||
          alevinOpts.protocol.end == bcEnd::THREE){
        // If the barcode sequence could be extracted, then this is set to true,
        // but the barcode sequence itself may still be invalid (e.g. contain `N` characters).
        // However, if extracted_barcode is false here, there is no hope to even recover the
        // barcode and we shouldn't attempt it.
        bool extracted_barcode = aut::extractBarcode(rp.first.seq, rp.second.seq, localProtocol, barcode);
        // If we could pull out something where the barcode sequence should have been
        // then continue to process it.
        if (extracted_barcode) {
          // if the barcode consisted of valid nucleotides, then seqOk is true
          // otherwise false
          seqOk =  aut::sequenceCheck(barcode, Sequence::BARCODE);
          if (not seqOk){
            // If the barcode contained invalid nucleotides
            // this attempts to replace the first one with an `A`.
            // If this returns true, there was only one `N` and we
            // replaced it; otherwise there was more than one `N`
            // and the barcode sequence should be treated as invalid.
            seqOk = aut::recoverBarcode(barcode);
          }
        }

        // If we have a barcode sequence, but not yet an index
        if (seqOk and (not barcodeIdx)) {
          // If we get here, we have a sequence-valid barcode.
          // Check if it is in the trBcs map.
          auto trIt = trBcs.find(barcode);

          // If it is, use that index
          if(trIt != trBcs.end()){
            barcodeIdx = trIt->second;
          } else{
            // If it's not, see if it's in the barcode map
            auto indIt = barcodeMap.find(barcode);
            // If so grab the representative and get its index
            if (indIt != barcodeMap.end()){
              barcode = indIt->second.front().first;
              auto trItLoc = trBcs.find(barcode);
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
          bool umi_ok = aut::extractUMI(rp.first.seq, rp.second.seq, localProtocol, umi);

          if ( !umi_ok ) {
            smallSeqs += 1;
          } else{
            alevin::types::AlevinUMIKmer umiIdx;
            bool isUmiIdxOk = umiIdx.fromChars(umi);

            if(isUmiIdxOk){
              jointHitGroup.setUMI(umiIdx.word(0));
	      if (writeQuasimappings) {
	      	extraBAMtags += "\tCB:Z:";
	      	extraBAMtags += barcode;
	      	extraBAMtags += "\tUR:Z:";
	      	extraBAMtags += umi;
	      }

              /*
              auto seq_len = rp.second.seq.size();
              if (alevinOpts.trimRight > 0) {
                if ( !tooShortRight ) {
                  //std::string sub_seq = rp.second.seq.substr(0, seq_len-alevinOpts.trimRight);
                  //auto rh = hitCollector(sub_seq, saSearcher, hcInfo);
                  readSubSeq = rp.second.seq.substr(0, seq_len-alevinOpts.trimRight);
                  auto rh = memCollector(readSubSeq, qc,
                                         true, // isLeft
                                         false // verbose
                                         );
                  //auto rh = hitCollector(readSubSeq, saSearcher, hcInfo);
                }
              } else {
              */
              readSubSeq = aut::getReadSequence(localProtocol, rp.first.seq, rp.second.seq, readBuffer);
              auto rh = tooShortRight ? false
                                      : memCollector(*readSubSeq, qc,
                                                     true, // isLeft
                                                     false // verbose
                                        );
              //}
              memCollector.findChains(*readSubSeq, hits,
                                      salmonOpts.fragLenDistMax,
                                      MateStatus::PAIRED_END_RIGHT,
                                      useChainingHeuristic, // heuristic chaining
                                      true, // isLeft
                                      false // verbose
                                      );

              pufferfish::util::joinReadsAndFilterSingle(hits, jointHits,
                                                         readSubSeq->length(),
                                                         memCollector.getConsensusFraction());
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
        if (tryAlign and !jointHits.empty()) {
          puffaligner.clear();
          puffaligner.getScoreStatus().reset();
          msi.clear(jointHits.size());

          size_t idx{0};
          bool is_multimapping = (jointHits.size() > 1);

          for (auto &&jointHit : jointHits) {
            // for alevin, currently, we need these to have a mate status of PAIRED_END_RIGHT
            jointHit.mateStatus = MateStatus::PAIRED_END_RIGHT;
            auto hitScore = puffaligner.calculateAlignments(*readSubSeq, jointHit, hctr, is_multimapping, salmonOpts.jointLog, false);
            bool validScore = (hitScore != invalidScore);
            numMappingsDropped += validScore ? 0 : 1;
            auto tid = qidx->getRefId(jointHit.tid);
            
            // NOTE: Here, we know that the read arising from the transcriptome is the "right"
            // read (read 2).  So we interpret compatibility in that context.
            // TODO: Make this code more generic and modular (account for the possibility of different library)
            // protocols or setups where the reads are not always "paired-end" and the transcriptomic read is not
            // always read 2 (@k3yavi).
            bool isCompat = (expectedLibraryFormat.strandedness == ReadStrandedness::U) or 
                            (jointHit.orphanClust()->isFw and (expectedLibraryFormat.strandedness == ReadStrandedness::AS)) or
                            (!jointHit.orphanClust()->isFw and (expectedLibraryFormat.strandedness == ReadStrandedness::SA));

            salmon::mapping_utils::updateRefMappings(tid, hitScore, isCompat, idx, transcripts, invalidScore, 
                                                     msi);
            ++idx;
          }

          bool bestHitDecoy = msi.haveOnlyDecoyMappings();
          if (msi.bestScore > invalidScore and !bestHitDecoy) {
            salmon::mapping_utils::filterAndCollectAlignments(jointHits,
                                                              readSubSeq->length(),
                                                              readSubSeq->length(),
                                                              false, // true for single-end false otherwise
                                                              tryAlign,
                                                              hardFilter,
                                                              salmonOpts.scoreExp,
                                                              salmonOpts.minAlnProb,
                                                              msi,
                                                              jointAlignments);
            if (!jointAlignments.empty()) {
              mapType = salmon::utils::MappingType::SINGLE_MAPPED;
            }
          } else {
            numDecoyFrags += bestHitDecoy ? 1 : 0;
            mapType = (bestHitDecoy) ? salmon::utils::MappingType::DECOY : salmon::utils::MappingType::UNMAPPED;
            if (bestHitDecoy) {
              salmon::mapping_utils::filterAndCollectAlignmentsDecoy(
                  jointHits, readSubSeq->length(),
                  readSubSeq->length(),
                  false, // true for single-end false otherwise
                  tryAlign, hardFilter, salmonOpts.scoreExp,
                  salmonOpts.minAlnProb, msi,
                  jointAlignments);
            } else {
              jointHitGroup.clearAlignments();
            }
          }
        } //end-if validate mapping

        if (writeQuasimappings) {
          writeAlignmentsToStream(rp, formatter, jointAlignments, sstream, true, true, extraBAMtags);
        }

        // We've kept decoy aignments around to this point so that we can
        // potentially write these alignments to the SAM file.  However, if 
        // we got to this point and only have decoy mappings, then clear the 
        // mappings here because none of the procesing below is relevant for 
        // decoys.
        if (mapType == salmon::utils::MappingType::DECOY) {
          jointHitGroup.clearAlignments();
        }

      if (writeUnmapped and mapType != salmon::utils::MappingType::SINGLE_MAPPED) {
        // If we have no mappings --- then there's nothing to do
        // unless we're outputting names for un-mapped reads
        unmappedNames << rp.first.name << ' ' << salmon::utils::str(mapType) << '\n';
      }

      validHits += jointAlignments.size();
      localNumAssignedFragments += (jointAlignments.size() > 0);
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
                                     readExp, fmCalc, rl, salmonOpts, hitLists,
                                     transcripts, clusterForest, fragLengthDist,
                                     numAssignedFragments, initialRound, burnedIn, maxZeroFrac);
  }

  if (maxZeroFrac > 5.0) {
    salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments",
                              maxZeroFrac);
  }

  mstats.numDecoyFragments += numDecoyFrags;
  readExp.updateShortFrags(shortFragStats);
}

template <typename AlnT, typename ProtocolT>
void sc_align_read_library(ReadExperimentT& readExp, 
                      ReadLibrary& rl, 
                      SalmonIndex* sidx, 
                      std::vector<Transcript>& transcripts, 
                      std::atomic<uint64_t>& numObservedFragments, // total number of reads we've looked at
                      std::atomic<uint64_t>& numAssignedFragments, // total number of assigned reads
                      std::atomic<uint32_t>& smallSeqs, 
                      std::atomic<uint32_t>& nSeqs, 
                      std::atomic<bool>& burnedIn, 
                      FragmentLengthDistribution& fragLengthDist, 
                      SalmonOpts& salmonOpts,
                      std::atomic<uint64_t>& num_chunks,
                      std::ofstream& rad_file,
                      std::ofstream& unmapped_bc_file,
                      std::mutex& fileMutex,
                      std::mutex& unmapped_bc_mutex,
                      std::mutex& ioMutex, 
                      size_t numThreads, 
                      std::vector<AlnGroupVec<AlnT>>& structureVec,
                      AlevinOpts<ProtocolT>& alevinOpts, 
                      MappingStatistics& mstats) {
  (void) burnedIn;
  std::vector<std::thread> threads;
  std::atomic<uint64_t> numValidHits{0};
  std::atomic<uint64_t> numUniqueMappings{0};
  rl.checkValid();
  auto indexType = sidx->indexType();
  std::unique_ptr<paired_parser> pairedParserPtr{nullptr};

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
    // Currently just a heuristic, a better way to do this would be
    // to have a joint thread-pool where threads could move between
    // being parsers and consumers.
    bool _did_modify = salmon::utils::configure_parsing(rl.mates1().size(), numThreads, numParsingThreads);
    //if (rl.mates1().size() > 1 and numThreads > 8) { numParsingThreads = 2; numThreads -= 1;}

    pairedParserPtr.reset(new paired_parser(rl.mates1(), rl.mates2(), numThreads, numParsingThreads, miniBatchSize));
    pairedParserPtr->start();

    auto processFunctor = [&](size_t i, auto* parserPtr, auto* index) {
     if(salmonOpts.qmFileName != "" and i == 0) {
       writeSAMHeader(*index, salmonOpts.qmLog);
     }
     auto threadFun = [&, i, parserPtr, index]() -> void {
       if (alevinOpts.sketch_mode) {
         process_reads_sc_sketch(
             parserPtr, readExp, rl, structureVec[i], numObservedFragments,
             numAssignedFragments, numUniqueMappings, numValidHits, smallSeqs, nSeqs, index,
             transcripts, fragLengthDist, salmonOpts, num_chunks, rad_file, unmapped_bc_file,
             fileMutex, unmapped_bc_mutex, ioMutex, alevinOpts, mstats
         );
       } else {
         process_reads_sc_align(
             parserPtr, readExp, rl, structureVec[i], numObservedFragments,
             numAssignedFragments, numValidHits, smallSeqs, nSeqs, index,
             transcripts, fragLengthDist, salmonOpts, num_chunks, rad_file, unmapped_bc_file,
             fileMutex, unmapped_bc_mutex, ioMutex, alevinOpts, mstats);
       }
     };
     threads.emplace_back(threadFun);
    };

    // True if we have a sparse index, false otherwise
    bool isSparse = sidx->isSparse();
    for (size_t i = 0; i < numThreads; ++i) {
      if (isSparse) {
        processFunctor(i, pairedParserPtr.get(), sidx->puffSparseIndex());
      } else {
        processFunctor(i, pairedParserPtr.get(), sidx->puffIndex());
      }
    } // End spawn all threads

    for (auto& t : threads) {
      t.join();
    }

    pairedParserPtr->stop();
  
    if (alevinOpts.sketch_mode) {
      salmonOpts.jointLog->info("Number uniquely mapped : {}", numUniqueMappings.load());
    }

    // DONE!
  }
}


template <typename AlnT, typename ProtocolT>
void processReadLibrary(
                        ReadExperimentT& readExp, 
                        ReadLibrary& rl, 
                        SalmonIndex* sidx,
                        std::vector<Transcript>& transcripts, 
                        ClusterForest& clusterForest,
                        std::atomic<uint64_t>& numObservedFragments, // total number of reads we've looked at
                        std::atomic<uint64_t>& numAssignedFragments, // total number of assigned reads
                        std::atomic<uint64_t>& upperBoundHits, // upper bound on # of mapped frags
                        std::atomic<uint32_t>& smallSeqs,
                        std::atomic<uint32_t>& nSeqs,
                        bool initialRound,
                        std::atomic<bool>& burnedIn, 
                        ForgettingMassCalculator& fmCalc,
                        FragmentLengthDistribution& fragLengthDist,
                        SalmonOpts& salmonOpts,
                        std::mutex& iomutex, 
                        size_t numThreads,
                        std::vector<AlnGroupVec<AlnT>>& structureVec,
                        AlevinOpts<ProtocolT>& alevinOpts,
                        SoftMapT& barcodeMap,
                        spp::sparse_hash_map<std::string, uint32_t>& trBcs, 
                        MappingStatistics& mstats) {

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
    auto processFunctor = [&](size_t i, auto* parserPtr, auto* index) {
     if(salmonOpts.qmFileName != "" and i == 0) {
       writeSAMHeader(*index, salmonOpts.qmLog);
     }
     auto threadFun = [&, i, parserPtr, index]() -> void {
                        processReadsQuasi(parserPtr,
                                          readExp, rl, structureVec[i],
                                          numObservedFragments, numAssignedFragments, numValidHits,
                                          upperBoundHits, smallSeqs, nSeqs, index, transcripts,
                                          fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                                          salmonOpts, iomutex, initialRound,
                                          burnedIn, alevinOpts, barcodeMap, trBcs,
                                          mstats);
     };
     threads.emplace_back(threadFun);
    };

    // True if we have a sparse index, false otherwise
    bool isSparse = sidx->isSparse();
    for (size_t i = 0; i < numThreads; ++i) {
      if (isSparse) {
        processFunctor(i, pairedParserPtr.get(), sidx->puffSparseIndex());
      } else {
        processFunctor(i, pairedParserPtr.get(), sidx->puffIndex());
      }
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
 * Selectively align the reads and write a PAM file out.
 */
template <typename AlnT, typename ProtocolT>
bool do_sc_align(ReadExperimentT& experiment,
                 SalmonOpts& salmonOpts,
                 MappingStatistics& mstats,
                 uint32_t numQuantThreads,
                 AlevinOpts<ProtocolT>& alevinOpts) {

  auto& refs = experiment.transcripts();
  size_t numTranscripts = refs.size();

  std::atomic<uint64_t> numObservedFragments{0};
  std::atomic<uint32_t> smallSeqs{0};
  std::atomic<uint32_t> nSeqs{0};

  auto jointLog = salmonOpts.jointLog;
  bool initialRound{true};
  uint32_t roundNum{0};

  std::mutex fileMutex;
  std::mutex unmapped_bc_mutex;
  std::mutex ioMutex;

  // RAD file path 
  boost::filesystem::path rad_file_path = salmonOpts.outputDirectory / "map.rad";
  boost::filesystem::path unmapped_bc_file_path = salmonOpts.outputDirectory / "unmapped_bc_count.bin";

  std::ofstream rad_file(rad_file_path.string());
  std::ofstream unmapped_bc_file(unmapped_bc_file_path.string());


  if (!rad_file.good()) {
    alevinOpts.jointLog->error("Could not open {} for writing.", rad_file_path.string());
    throw std::runtime_error("error creating output file.");
  }

  if (!unmapped_bc_file.good()) {
    alevinOpts.jointLog->error("Could not open {} for writing.", unmapped_bc_file_path.string());
    throw std::runtime_error("error creating output file.");
  }



  BasicBinWriter bw;
  //  RADHeader
  RADHeader rh;
  // right now, all formats are effectively single-end
  rh.is_paired(false);
  for (auto& r : refs) {
    rh.add_refname(r.RefName);
  }
  rh.dump_to_bin(bw);
  
  // where we will write the number of chunks when we know
  // how many there are
  auto chunk_offset = bw.num_bytes() - sizeof(uint64_t);
  std::atomic<uint64_t> num_chunks{0};

  // ### start of tags 

  // Tags we will have
  // write the tag meta-information section

  // File-level tag description 
  uint16_t file_level_tags{2};
  bw << file_level_tags;

  // cblen
  uint8_t type_id{2};
  bw << std::string("cblen");
  bw << type_id;

  bw << std::string("ulen");
  bw << type_id;

  // read-level tag description
  uint16_t read_level_tags{2};
  bw << read_level_tags;

  // barcode
  bw << std::string("b");
  if ( alevinOpts.protocol.barcodeLength > 32 ) {
    type_id = 8;
  } else if ( alevinOpts.protocol.barcodeLength > 16 ) {
    // 17 - 32 bases
    type_id = 4;
  } else {
    // <= 16 bases 
    type_id = 3;
  }
  bw << type_id;

  // umi
  bw << std::string("u");
  if ( alevinOpts.protocol.umiLength > 32 ) {
    type_id = 8;
  } else if ( alevinOpts.protocol.umiLength > 16 ) {
    // 17 - 32 bases
    type_id = 4;
  } else {
    // <= 16 bases 
    type_id = 3;
  }
  bw << type_id;

  // alignment-level tag description
  uint16_t aln_level_tags{1};
  bw << aln_level_tags;
  // we maintain orientation
  //bw << std::string("orientation");
  //type_id = 1;
  //bw << type_id;

  // and reference id
  bw << std::string("compressed_ori_refid");
  type_id = 3;
  bw << type_id;

  // ### end of tag definitions

  // the actual file-level tags
  bw << static_cast<uint16_t>(alevinOpts.protocol.barcodeLength);
  bw << static_cast<uint16_t>(alevinOpts.protocol.umiLength);

  rad_file << bw;
  bw.clear();

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

  auto processReadLibraryCallback =
    [&](ReadLibrary& rl, SalmonIndex* sidx,
        std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
        FragmentLengthDistribution& fragLengthDist,
        std::atomic<uint64_t>& numAssignedFragments, size_t numQuantThreads,
        std::atomic<bool>& burnedIn) -> void {
    (void) clusterForest;
    sc_align_read_library<AlnT>(experiment, rl, sidx, transcripts, 
                             numObservedFragments, numAssignedFragments,
                             smallSeqs, nSeqs, burnedIn, 
                             fragLengthDist, salmonOpts,
                             num_chunks,
                             rad_file,
                             unmapped_bc_file,
                             fileMutex,
                             unmapped_bc_mutex,
                             ioMutex, numQuantThreads, groupVec,
                             alevinOpts, mstats);
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
      auto* sidx = experiment.getIndex();
      bool isSparse = sidx->isSparse();
      size_t minK = (isSparse) ? sidx->puffSparseIndex()->k() : sidx->puffIndex()->k();
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
  salmonOpts.jointLog->info("Selectively-aligned {} total fragments out of {}", experiment.numAssignedFragments(), numObservedFragments.load() );
  salmonOpts.jointLog->info("Number of fragments discarded because they are best-mapped to decoys : {:n}",
                            mstats.numDecoyFragments.load());

  if (smallSeqs > 100) {
    jointLog->warn("Found {} reads with CB+UMI length smaller than expected.\n"
                   "Please report on github if this number is too large", smallSeqs);
  }
  if (nSeqs > 100) {
    jointLog->warn("Found {} reads with `N` in the UMI sequence and ignored the reads.\n"
                   "Please report on github if this number is too large", nSeqs);
  }

  rad_file.seekp(chunk_offset);
  uint64_t nc = num_chunks.load();
  rad_file.write(reinterpret_cast<char*>(&nc), sizeof(nc));

  rad_file.close();

  // We want to check if the RAD file stream was written to properly.
  // While we likely would have caught this earlier, it is possible the 
  // badbit may not be set until the stream actually flushes (perhaps even 
  // at close), so we check here one final time that the status of the 
  // stream is as expected.
  // see : https://stackoverflow.com/questions/28342660/error-handling-in-stdofstream-while-writing-data
  if ( !rad_file ) {
    jointLog->critical("The RAD file stream had an invalid status after close; so some operation(s) may"
                       "have failed!\nA common cause for this is lack of output disk space.\n"
                       "Consequently, the output may be corrupted.");
    jointLog->flush();
    return false;
  }

  jointLog->info("finished sc_align()");
  jointLog->flush();
  return true;
}
                     
/**
 *  Quantify the targets given in the file `transcriptFile` using the
 *  reads in the given set of `readLibraries`, and write the results
 *  to the file `outputFile`.  The reads are assumed to be in the format
 *  specified by `libFmt`.
 *
 */
template <typename AlnT, typename ProtocolT>
void quantifyLibrary(ReadExperimentT& experiment,
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
                             ioMutex,
                             numQuantThreads, groupVec,
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
      auto* sidx = experiment.getIndex();
      bool isSparse = sidx->isSparse();
      size_t minK = (isSparse) ? sidx->puffSparseIndex()->k() : sidx->puffIndex()->k();
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
int alevin_sc_align(AlevinOpts<ProtocolT>& aopt,
                    SalmonOpts& sopt,
                    boost::program_options::parsed_options& orderedOptions,
                    std::unique_ptr<SalmonIndex>& salmonIndex){
  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  try {
    auto jointLog = aopt.jointLog;
    auto indexDirectory = sopt.indexDirectory;
    auto outputDirectory = sopt.outputDirectory;

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

    if(!salmonIndex)
      salmonIndex = checkLoadIndex(indexDirectory, sopt.jointLog);
    auto idxType = salmonIndex->indexType();

    MappingStatistics mstats;
    ReadExperimentT experiment(readLibraries, salmonIndex.get(), sopt);

    // We currently do not support decoy sequence in the 
    // --justAlign or --sketch modes, so check that the 
    // index is free of decoys, otherwise complain and exit.
    bool index_has_decoys = false;
    auto* idx = experiment.getIndex();
    // if this is a sparse index
    if (idx->isSparse()) {
      auto* pfi = idx->puffSparseIndex();
      auto fdei = pfi->firstDecoyEncodedIndex();
      index_has_decoys = (fdei != std::numeric_limits<decltype(fdei)>::max());
    } else {
      // otherwise this is a dense index
      auto* pfi = idx->puffIndex();
      auto fdei = pfi->firstDecoyEncodedIndex();
      index_has_decoys = (fdei != std::numeric_limits<decltype(fdei)>::max());
    }

    if (index_has_decoys) {
      jointLog->error("The provided index has decoy sequences.  However, decoy "
                      "sequences are not currently supported in the \"rad\" or "
                      "\"sketch\" modes.  Please provide an index without decoy "
                      "sequence to continue.");
      jointLog->flush();
      return 1;
    }

    // This will be the class in charge of maintaining our
    // rich equivalence classes
    experiment.equivalenceClassBuilder().setMaxResizeThreads(sopt.maxHashResizeThreads);
    experiment.equivalenceClassBuilder().start();

    auto indexType = experiment.getIndex()->indexType();

    sopt.allowOrphans = true;
    sopt.useQuasi = true;
    
    // lose one thread for parsing
    // TODO: is it reasonable to assume 
    // that the work done by the parsing thread
    // will be substantial enough to count if 
    // the remaining mapping thread count is 
    // too small?
    if(sopt.numThreads > 1){
      sopt.numThreads -= 1;
    }

    if (aopt.protocol.barcodeLength <= 32) {
      alevin::types::AlevinCellBarcodeKmer::k(aopt.protocol.barcodeLength);
    }

    bool mapping_ok = do_sc_align<QuasiAlignment>(experiment, sopt,
                                                  mstats, sopt.numThreads, aopt);

    // write meta-information about the run
    GZipWriter gzw(outputDirectory, jointLog);
    sopt.runStopTime = salmon::utils::getCurrentTimeAsString();
    gzw.writeMetaFryMode(sopt, experiment, mstats);

    if ( !mapping_ok ) {
      std::cerr << "Error : There was an error during mapping; please check the log for further details.\n";
      std::exit(1);
    }
 } catch (po::error& e) {
    std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (const spdlog::spdlog_ex& ex) {
    std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
    std::exit(1);
  } catch (std::exception& e) {
    std::cerr << "Exception : [" << e.what() << "]\n";
    std::cerr << " alevin (sc-align) was invoked improperly.\n";
    std::cerr << "For usage information, try "
              << " alevin --help\nExiting.\n";
    std::exit(1);
  }
  return 0;
}

template <typename ProtocolT>
int alevinQuant(AlevinOpts<ProtocolT>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
                std::unique_ptr<SalmonIndex>& salmonIndex){
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
    if(!salmonIndex)
      salmonIndex = checkLoadIndex(indexDirectory, sopt.jointLog);
    auto idxType = salmonIndex->indexType();

    MappingStatistics mstats;
    ReadExperimentT experiment(readLibraries, salmonIndex.get(), sopt);
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
    quantifyLibrary<QuasiAlignment>(experiment, sopt,
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
int alevin_sc_align(AlevinOpts<apt::DropSeq>& aopt, SalmonOpts& sopt,
                    boost::program_options::parsed_options& orderedOptions,
                    std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevinQuant(AlevinOpts<apt::DropSeq>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template
int alevin_sc_align(AlevinOpts<apt::CITESeq>& aopt, SalmonOpts& sopt,
                    boost::program_options::parsed_options& orderedOptions,
                    std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::CITESeq>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::InDropV2>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::InDropV2>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::ChromiumV3>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::ChromiumV3>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::Chromium>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::Chromium>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::Gemcode>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::Gemcode>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::CELSeq>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::CELSeq>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::CELSeq2>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::CELSeq2>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::SplitSeqV1>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::SplitSeqV1>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::SplitSeqV2>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);
template int
alevinQuant(AlevinOpts<apt::SplitSeqV2>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::QuartzSeq2>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevinQuant(AlevinOpts<apt::QuartzSeq2>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template 
int alevin_sc_align(AlevinOpts<apt::SciSeq3>& aopt,
                    SalmonOpts& sopt,
                    boost::program_options::parsed_options& orderedOptions,
                    std::unique_ptr<SalmonIndex>& salmonIndex);
template
int alevinQuant(AlevinOpts<apt::SciSeq3>& aopt,
                SalmonOpts& sopt,
                SoftMapT& barcodeMap,
                TrueBcsT& trueBarcodes,
                spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
                spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
                boost::program_options::parsed_options& orderedOptions,
                CFreqMapT& freqCounter,
                size_t numLowConfidentBarcode,
                std::unique_ptr<SalmonIndex>& salmonIndex);


template int
alevin_sc_align(AlevinOpts<apt::Custom>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevinQuant(AlevinOpts<apt::Custom>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevin_sc_align(AlevinOpts<apt::CustomGeometry>& aopt, SalmonOpts& sopt,
                boost::program_options::parsed_options& orderedOptions,
                std::unique_ptr<SalmonIndex>& salmonIndex);

template int
alevinQuant(AlevinOpts<apt::CustomGeometry>& aopt, SalmonOpts& sopt,
            SoftMapT& barcodeMap, TrueBcsT& trueBarcodes,
            spp::sparse_hash_map<uint32_t, uint32_t>& txpToGeneMap,
            spp::sparse_hash_map<std::string, uint32_t>& geneIdxMap,
            boost::program_options::parsed_options& orderedOptions,
            CFreqMapT& freqCounter, size_t numLowConfidentBarcode,
            std::unique_ptr<SalmonIndex>& salmonIndex);
