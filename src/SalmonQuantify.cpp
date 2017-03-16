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

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "ksort.h"
#include "kvec.h"
#include "utils.h"
}

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
//#include "TextBootstrapWriter.hpp"

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

template <typename AlnT> using AlnGroupVec = std::vector<AlignmentGroup<AlnT>>;

template <typename AlnT>
using AlnGroupVecRange =
    boost::iterator_range<typename AlnGroupVec<AlnT>::iterator>;

#define __MOODYCAMEL__
#if defined(__MOODYCAMEL__)
template <typename AlnT>
using AlnGroupQueue = moodycamel::ConcurrentQueue<AlignmentGroup<AlnT>*>;
#else
template <typename AlnT>
using AlnGroupQueue = tbb::concurrent_queue<AlignmentGroup<AlnT>*>;
#endif

#include "LightweightAlignmentDefs.hpp"

template <typename AlnT>
void processMiniBatch(ReadExperiment& readExp, ForgettingMassCalculator& fmCalc,
                      uint64_t firstTimestepOfRound, ReadLibrary& readLib,
                      const SalmonOpts& salmonOpts,
                      AlnGroupVecRange<AlnT> batchHits,
                      std::vector<Transcript>& transcripts,
                      ClusterForest& clusterForest,
                      FragmentLengthDistribution& fragLengthDist,
                      BiasParams& observedBiasParams,
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
  EquivalenceClassBuilder& eqBuilder = readExp.equivalenceClassBuilder();

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
                // Threading is hard!  It's possible that an update to the PMF snuck in between when we asked to cache the CMF and when the
                // "burnedIn" variable was last seen as false.
                log->info("reference length = {}, CMF[refLen] = {}, fragLen = {}, PMF[fragLen] = {}", refLength, std::exp(refLengthCM), aln.fragLength(), std::exp(lenProb));
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

	  /*
	  double logAlignCompatProb =
	      (useReadCompat) ? (salmon::utils::logAlignFormatProb(
				    aln.libFormat(), expectedLibraryFormat,
				    static_cast<int32_t>(aln.pos), aln.fwd,
                                    aln.mateStatus, salmonOpts.incompatPrior))
                              : LOG_1;
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
          // DEC 9
          if (aln.mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED and !noLengthCorrection) {
            startPosProb = (flen <= refLength) ? -std::log(refLength - flen + 1) : salmon::math::LOG_EPSILON;
          }

          double fragStartLogNumerator{salmon::math::LOG_1};
          double fragStartLogDenominator{salmon::math::LOG_1};

          auto hitPos = aln.hitPos();
          /** NOTE: no more FSPD
          if (useFSPD and burnedIn and hitPos < refLength) {
            auto& fragStartDist = fragStartDists[transcript.lengthClassIndex()];
            // Get the log(numerator) and log(denominator) for the fragment
            // start position
            // probability.
            bool nonZeroProb = fragStartDist.logNumDenomMass(
                hitPos, refLength, logRefLength, fragStartLogNumerator,
                fragStartLogDenominator);
            // Set the overall probability.
            startPosProb = (nonZeroProb)
                               ? fragStartLogNumerator - fragStartLogDenominator
                               : salmon::math::LOG_0;
          }
          **/

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
          if (p >= transcript.RefLength) {
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
              posFW = posFW >= transcript.RefLength ? transcript.RefLength - 1
                                                    : posFW;
              posRC = posRC < 0 ? 0 : posRC;
              posRC = posRC >= transcript.RefLength ? transcript.RefLength - 1
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
            pos = pos >= transcript.RefLength ? transcript.RefLength - 1 : pos;
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
                if (start >= 0 and stop < transcript.RefLength) {
                    bool valid{false};
                    auto desc = transcript.gcDesc(start, stop, valid);
                    if (valid) { observedGCMass.inc(desc, aln.logProb); }
                }
            } else if(expectedLibraryFormat.type == ReadType::SINGLE_END) { 
	      // Both expected and observed should be single end here
                // For single-end reads, simply assume that every fragment
                // has a length equal to the conditional mean (given the 
                // current transcript's length).
                auto cmeans = readExp.condMeans();
                auto cmean = static_cast<int32_t>((transcript.RefLength >= cmeans.size()) ? cmeans.back() : cmeans[transcript.RefLength]);
                int32_t start = aln.fwd ? aln.pos : std::max(0, aln.pos - cmean);
                int32_t stop = start + cmean;
                // WITH CONTEXT
                if (start >= 0 and stop < transcript.RefLength) {
                  bool valid{false};
                  auto desc = transcript.gcDesc(start, stop, valid);
                  if (valid) {observedGCMass.inc(desc, aln.logProb);}
                }
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
template <typename RapMapIndexT>
void processReadsQuasi(
    paired_parser* parser, ReadExperiment& readExp, ReadLibrary& rl,
    AlnGroupVec<SMEMAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    RapMapIndexT* idx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
    mem_opt_t* memOptions, SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache) {

  // ERROR
  salmonOpts.jointLog->error("MEM-mapping cannot be used with the Quasi index "
                             "--- please report this bug on GitHub");
  std::exit(1);
}

template <typename RapMapIndexT>
void processReadsQuasi(
    single_parser* parser, ReadExperiment& readExp, ReadLibrary& rl,
    AlnGroupVec<SMEMAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    RapMapIndexT* sidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
    mem_opt_t* memOptions, SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache) {
  // ERROR
  salmonOpts.jointLog->error("MEM-mapping cannot be used with the Quasi index "
                             "--- please report this bug on GitHub");
  std::exit(1);
}

template <typename RapMapIndexT>
void processReadsQuasi(
    paired_parser* parser, ReadExperiment& readExp, ReadLibrary& rl,
    AlnGroupVec<QuasiAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    RapMapIndexT* qidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
    mem_opt_t* memOptions, SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache) {
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

  auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
  size_t minK = rapmap::utils::my_mer::k();

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};
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
      localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      jointHitGroup.clearAlignments();
      auto& jointHits = jointHitGroup.alignments();
      leftHits.clear();
      rightHits.clear();
      mapType = salmon::utils::MappingType::UNMAPPED;

      bool lh = tooShortLeft ? false : hitCollector(rp.first.seq,
                                                    leftHits, saSearcher,
                                                    MateStatus::PAIRED_END_LEFT,
                                                    consistentHits);

      bool rh = tooShortRight ? false : hitCollector(rp.second.seq,
                                   rightHits, saSearcher,
                                   MateStatus::PAIRED_END_RIGHT, 
                                   consistentHits);

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
                (startPos1 > 0 and startPos1 < t.RefLength) and
                (startPos2 > 0 and startPos2 < t.RefLength)) {

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
                       t.RefLength) and
                  (startPos2 >= readBias2.contextBefore(read2RC) and
                   startPos2 + readBias2.contextAfter(read2RC) < t.RefLength)) {

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
            h.format = salmon::utils::hitType(h.pos, !h.fwd);
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
          fmt::print(stderr, "\033[A\r\r{}processed{} {} {}fragments{}\n",
                     green, red, numObservedFragments, green, RESET_COLOR);
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
        numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac);
  }

  if (maxZeroFrac > 0.0) {
      salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments", 
                                maxZeroFrac);
  }

  readExp.updateShortFrags(shortFragStats);
}

// SINGLE END

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename RapMapIndexT>
void processReadsQuasi(
    single_parser* parser, ReadExperiment& readExp, ReadLibrary& rl,
    AlnGroupVec<QuasiAlignment>& structureVec,
    std::atomic<uint64_t>& numObservedFragments,
    std::atomic<uint64_t>& numAssignedFragments,
    std::atomic<uint64_t>& validHits, std::atomic<uint64_t>& upperBoundHits,
    RapMapIndexT* qidx, std::vector<Transcript>& transcripts,
    ForgettingMassCalculator& fmCalc, ClusterForest& clusterForest,
    FragmentLengthDistribution& fragLengthDist, BiasParams& observedBiasParams,
    mem_opt_t* memOptions, SalmonOpts& salmonOpts, double coverageThresh,
    std::mutex& iomutex, bool initialRound, std::atomic<bool>& burnedIn,
    volatile bool& writeToCache) {
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
  spdlog::logger* unmappedLogger = (writeUnmapped) ? salmonOpts.unmappedLog.get() : nullptr;

  auto& readBiasFW = observedBiasParams.seqBiasModelFW;
  auto& readBiasRC = observedBiasParams.seqBiasModelRC;
  Mer context;

  const char* txomeStr = qidx->seq.c_str();

  auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
  size_t minK = rapmap::utils::my_mer::k();

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};

  bool tooManyHits{false};
  size_t readLen{0};
  size_t maxNumHits{salmonOpts.maxReadOccs};
  bool consistentHits{salmonOpts.consistentHits};
  bool quiet{salmonOpts.quiet};

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
  rapmap::utils::HitCounters hctr;
  
  SingleAlignmentFormatter<RapMapIndexT*> formatter(qidx);
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
      readLen = rp.seq.length();
      tooShort = (readLen < minK);
      tooManyHits = false;
      localUpperBoundHits = 0;
      auto& jointHitGroup = structureVec[i];
      auto& jointHits = jointHitGroup.alignments();
      jointHitGroup.clearAlignments();

      bool lh =
          tooShort ? false
          : hitCollector(rp.seq,
                                  jointHits, saSearcher,
                                  MateStatus::SINGLE_END, consistentHits);

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
          if (startPos > 0 and startPos < t.RefLength) {
            auto& readBias = (h.fwd) ? readBiasFW : readBiasRC;
            const char* txpStart = t.Sequence();
            const char* readStart = txpStart + startPos;
            const char* txpEnd = txpStart + t.RefLength;

            bool success{false};
            // If the context exists around the read, add it to the observed
            // read start sequences.
            if (startPos >= readBias.contextBefore(!h.fwd) and
                startPos + readBias.contextAfter(!h.fwd) < t.RefLength) {
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
        }
      }
      
      if (writeQuasimappings) {
          rapmap::utils::writeAlignmentsToStream(rp, formatter,
                                                 hctr, jointHits, sstream);
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
          fmt::print(stderr, "\033[A\r\r{}processed{} {} {}fragments{}\n",
                     green, red, numObservedFragments, green, RESET_COLOR);
          fmt::print(stderr, "hits: {}; hits per frag:  {}", validHits,
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
    AlnGroupVecRange<QuasiAlignment> hitLists = boost::make_iterator_range(
        structureVec.begin(), structureVec.begin() + rangeSize);
    processMiniBatch<QuasiAlignment>(
        readExp, fmCalc, firstTimestepOfRound, rl, salmonOpts, hitLists,
        transcripts, clusterForest, fragLengthDist, observedBiasParams,
        numAssignedFragments, eng, initialRound, burnedIn, maxZeroFrac);
  }
  readExp.updateShortFrags(shortFragStats);

  if (maxZeroFrac > 0.0) {
      salmonOpts.jointLog->info("Thread saw mini-batch with a maximum of {0:.2f}\% zero probability fragments", 
                                maxZeroFrac);
  }
}

/// DONE QUASI

template <typename AlnT>
void processReadLibrary(
    ReadExperiment& readExp, ReadLibrary& rl, SalmonIndex* sidx,
    std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
    std::atomic<uint64_t>&
        numObservedFragments, // total number of reads we've looked at
    std::atomic<uint64_t>&
        numAssignedFragments,              // total number of assigned reads
    std::atomic<uint64_t>& upperBoundHits, // upper bound on # of mapped frags
    bool initialRound,
    std::atomic<bool>& burnedIn, ForgettingMassCalculator& fmCalc,
    FragmentLengthDistribution& fragLengthDist, mem_opt_t* memOptions,
    SalmonOpts& salmonOpts, double coverageThresh, bool greedyChain,
    std::mutex& iomutex, size_t numThreads,
    std::vector<AlnGroupVec<AlnT>>& structureVec, volatile bool& writeToCache) {

  std::vector<std::thread> threads;

  std::atomic<uint64_t> numValidHits{0};
  rl.checkValid();

  auto indexType = sidx->indexType();

  std::unique_ptr<paired_parser> pairedParserPtr{nullptr};
  std::unique_ptr<single_parser> singleParserPtr{nullptr};

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
    if (rl.mates1().size() > 1 and numThreads > 8) { numParsingThreads = 2; }
    pairedParserPtr.reset(new paired_parser(rl.mates1(), rl.mates2(), numThreads, numParsingThreads, miniBatchSize));
    pairedParserPtr->start();
    
    switch (indexType) {
    case SalmonIndexType::FMD: {
      for (int i = 0; i < numThreads; ++i) {
        // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
        // change value before the lambda below is evaluated --- crazy!
        auto threadFun = [&, i]() -> void {
          processReadsMEM<paired_parser, TranscriptHitList>(
              pairedParserPtr.get(), readExp, rl, structureVec[i],
              numObservedFragments, numAssignedFragments, numValidHits,
              upperBoundHits, sidx, transcripts, fmCalc, clusterForest,
              fragLengthDist, observedBiasParams[i], memOptions, salmonOpts,
              coverageThresh, iomutex, initialRound, burnedIn, writeToCache);
        };
        threads.emplace_back(threadFun);
      }
      break;
    case SalmonIndexType::QUASI: {
      // True if we have a 64-bit SA index, false otherwise
      bool largeIndex = sidx->is64BitQuasi();
      bool perfectHashIndex = sidx->isPerfectHashQuasi();
      for (int i = 0; i < numThreads; ++i) {
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
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
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
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
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
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
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
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
            };
            threads.emplace_back(threadFun);
          }

        } // End spawn current thread

      } // End spawn all threads
    }   // End Quasi index
    break;
    } // end switch
    }
    for (int i = 0; i < numThreads; ++i) {
      threads[i].join();
    }

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

  } // ------ Single-end --------
  else if (rl.format().type == ReadType::SINGLE_END) {

    uint32_t numParsingThreads{1};
    // HACK!
    if (rl.unmated().size() > 1 and numThreads > 8) { numParsingThreads = 2; }
    singleParserPtr.reset(new single_parser(rl.unmated(), numThreads, numParsingThreads, miniBatchSize));
    singleParserPtr->start();
    switch (indexType) {
    case SalmonIndexType::FMD: {
      for (int i = 0; i < numThreads; ++i) {
        // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
        // change value before the lambda below is evaluated --- crazy!
        auto threadFun = [&, i]() -> void {
          processReadsMEM<single_parser, TranscriptHitList>(
              singleParserPtr.get(), readExp, rl, structureVec[i],
              numObservedFragments, numAssignedFragments, numValidHits,
              upperBoundHits, sidx, transcripts, fmCalc, clusterForest,
              fragLengthDist, observedBiasParams[i], memOptions, salmonOpts,
              coverageThresh, iomutex, initialRound, burnedIn, writeToCache);
        };
        threads.emplace_back(threadFun);
      }
    } break;

    case SalmonIndexType::QUASI: {
      // True if we have a 64-bit SA index, false otherwise
      bool largeIndex = sidx->is64BitQuasi();
      bool perfectHashIndex = sidx->isPerfectHashQuasi();
      for (int i = 0; i < numThreads; ++i) {

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
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
            };
            threads.emplace_back(threadFun);
          } else { // Dense Hash
              if (salmonOpts.qmFileName != "" and i == 0) {
                  rapmap::utils::writeSAMHeader(*(sidx->quasiIndex64()), salmonOpts.qmLog);
              }

            auto threadFun = [&, i]() -> void {
              processReadsQuasi<RapMapSAIndex<int64_t, DenseHash<int64_t>>>(
                  singleParserPtr.get(), readExp, rl, structureVec[i],
                  numObservedFragments, numAssignedFragments, numValidHits,
                  upperBoundHits, sidx->quasiIndex64(), transcripts, fmCalc,
                  clusterForest, fragLengthDist, observedBiasParams[i],
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
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
                  singleParserPtr.get(), readExp, rl, structureVec[i],
                  numObservedFragments, numAssignedFragments, numValidHits,
                  upperBoundHits, sidx->quasiIndexPerfectHash32(), transcripts,
                  fmCalc, clusterForest, fragLengthDist, observedBiasParams[i],
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
            };
            threads.emplace_back(threadFun);
          } else { // Dense Hash
              if (salmonOpts.qmFileName != "" and i == 0) {
                  rapmap::utils::writeSAMHeader(*(sidx->quasiIndex32()), salmonOpts.qmLog);
              }

            auto threadFun = [&, i]() -> void {
              processReadsQuasi<RapMapSAIndex<int32_t, DenseHash<int32_t>>>(
                  singleParserPtr.get(), readExp, rl, structureVec[i],
                  numObservedFragments, numAssignedFragments, numValidHits,
                  upperBoundHits, sidx->quasiIndex32(), transcripts, fmCalc,
                  clusterForest, fragLengthDist, observedBiasParams[i],
                  memOptions, salmonOpts, coverageThresh, iomutex, initialRound,
                  burnedIn, writeToCache);
            };
            threads.emplace_back(threadFun);
          }

        } // End spawn current thread

      } // End spawn all threads
    }   // End Quasi index
    break;
    }
    for (int i = 0; i < numThreads; ++i) {
      threads[i].join();
    }

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

    /* OLD SINGLE END BIAS
    // Set the global distribution based on the sum of local
    // distributions.
    for (auto& gcp : observedBiasParams) {
      auto& fw = readExp.readBiasModelObserved(salmon::utils::Direction::FORWARD);
      auto& rc =
          readExp.readBiasModelObserved(salmon::utils::Direction::REVERSE_COMPLEMENT);

      auto& fwloc = gcp.seqBiasModelFW;
      auto& rcloc = gcp.seqBiasModelRC;
      fw.combineCounts(fwloc);
      rc.combineCounts(rcloc);

      // positional biases
      auto& posBiasesFW = readExp.posBias(salmon::utils::Direction::FORWARD);
      auto& posBiasesRC =
          readExp.posBias(salmon::utils::Direction::REVERSE_COMPLEMENT);
      for (size_t i = 0; i < posBiasesFW.size(); ++i) {
        posBiasesFW[i].combine(gcp.posBiasFW[i]);
        posBiasesRC[i].combine(gcp.posBiasRC[i]);
      }
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
    END OLD SINGLE-END BIAS */


  } // ------ END Single-end --------
}

/**
  *  Quantify the targets given in the file `transcriptFile` using the
  *  reads in the given set of `readLibraries`, and write the results
  *  to the file `outputFile`.  The reads are assumed to be in the format
  *  specified by `libFmt`.
  *
  */
template <typename AlnT>
void quantifyLibrary(ReadExperiment& experiment, bool greedyChain,
                     mem_opt_t* memOptions, SalmonOpts& salmonOpts,
                     double coverageThresh, uint32_t numQuantThreads) {

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
                               fragLengthDist, memOptions, salmonOpts,
                               coverageThresh, greedyChain, ioMutex,
                               numQuantThreads, groupVec, writeToCache);

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

  jointLog->info("Mapping rate = {}\%\n",
                 experiment.effectiveMappingRate() * 100.0);
  jointLog->info("finished quantifyLibrary()");
}

int salmonQuantify(int argc, char* argv[]) {
  using std::cerr;
  using std::vector;
  using std::string;
  namespace bfs = boost::filesystem;
  namespace po = boost::program_options;

  bool optChain{false};
  int32_t numBiasSamples{0};

  SalmonOpts sopt;
  mem_opt_t* memOptions = mem_opt_init();
  memOptions->split_factor = 1.5;

  sopt.numThreads = std::thread::hardware_concurrency();

  double coverageThresh;
  vector<string> unmatedReadFiles;
  vector<string> mate1ReadFiles;
  vector<string> mate2ReadFiles;

  po::options_description generic("\n"
                                  "basic options");
  generic.add_options()("version,v", "print version string")
    (
      "help,h", "produce help message")
    (
      "index,i", po::value<string>()->required(),
      "Salmon index")("libType,l", po::value<std::string>()->required(),
                      "Format string describing the library type")
    (
      "unmatedReads,r",
      po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
      "List of files containing unmated reads of (e.g. single-end reads)")(
      "mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
      "File containing the #1 mates")
    (
      "mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
      "File containing the #2 mates")
    (

      "output,o", po::value<std::string>()->required(),
      "Output quantification file.")
    (
      "allowOrphans",
      po::bool_switch(&(sopt.allowOrphans))->default_value(false),
      "Consider orphaned reads as valid hits when "
      "performing lightweight-alignment.  This option will increase "
      "sensitivity (allow more reads to map and "
      "more transcripts to be detected), but may decrease specificity as "
      "orphaned alignments are more likely "
      "to be spurious -- this option is *always* set to true when using "
      "quasi-mapping.")
    (
     "seqBias",
     po::bool_switch(&(sopt.biasCorrect))->default_value(false),
     "Perform sequence-specific bias correction.")
    (
      "gcBias", po::bool_switch(&(sopt.gcBiasCorrect))->default_value(false),
      "[beta] Perform fragment GC bias correction")
    (
      "threads,p",
      po::value<uint32_t>(&(sopt.numThreads))->default_value(sopt.numThreads),
      "The number of threads to use concurrently.")
    (
      "incompatPrior",
      po::value<double>(&(sopt.incompatPrior))->default_value(1e-20),
      "This option "
      "sets the prior probability that an alignment that disagrees with the "
      "specified "
      "library type (--libType) results from the true fragment origin.  "
      "Setting this to 0 "
      "specifies that alignments that disagree with the library type should be "
      "\"impossible\", "
      "while setting it to 1 says that alignments that disagree with the "
      "library type are no "
      "less likely than those that do")
    (
      "geneMap,g", po::value<string>(),
      "File containing a mapping of transcripts to genes.  If this file is "
      "provided "
      "Salmon will output both quant.sf and quant.genes.sf files, where the "
      "latter "
      "contains aggregated gene-level abundance estimates.  The transcript to "
      "gene mapping "
      "should be provided as either a GTF file, or a in a simple tab-delimited "
      "format "
      "where each line contains the name of a transcript and the gene to which "
      "it belongs "
      "separated by a tab.  The extension of the file is used to determine how "
      "the file "
      "should be parsed.  Files ending in \'.gtf\', \'.gff\' or \'.gff3\' are assumed to "
      "be in GTF "
      "format; files with any other extension are assumed to be in the simple "
      "format. In GTF / GFF format, the \"transcript_id\" is assumed to contain the "
      "transcript identifier and the \"gene_id\" is assumed to contain the corresponding "
      "gene identifier.")
  (
   "writeMappings,z", po::value<string>(&sopt.qmFileName)->default_value("")->implicit_value("-"),
   "If this option is provided, then the quasi-mapping results will be written out in SAM-compatible "
   "format.  By default, output will be directed to stdout, but an alternative file name can be "
   "provided instead.")
  (
   "meta", po::bool_switch(&(sopt.meta))->default_value(false),
   "If you're using Salmon on a metagenomic dataset, consider setting this flag to disable parts of the "
   "abundance estimation model that make less sense for metagenomic data."
  );

  sopt.noRichEqClasses = false;
  // mapping cache has been deprecated
  sopt.disableMappingCache = true;

  po::options_description advanced("\n"
                                   "advanced options");
  advanced.add_options()
      /*
      ("disableMappingCache",
      po::bool_switch(&(sopt.disableMappingCache))->default_value(false),
      "Setting this option disables the creation and use "
                                          "of the \"mapping cache\" file.  The
      mapping cache can speed up quantification significantly for smaller read "
                                          "libraries (i.e. where the number of
      mapped fragments is less than the required number of observations).
      However, "
                                          "for very large read libraries, the
      mapping cache is unnecessary, and disabling it may allow salmon to more
      effectively "
                                          "make use of a very large number of
      threads.")
      */
    ("alternativeInitMode", po::bool_switch(&(sopt.alternativeInitMode))->default_value(false),
     "[Experimental]: Use an alternative strategy (rather than simple interpolation between) the "
     "online and uniform abundance estimates to initalize the EM / VBEM algorithm."
     )
    (
     "auxDir", po::value<std::string>(&(sopt.auxDir))->default_value("aux_info"),
     "The sub-directory of the quantification directory where auxiliary "
     "information "
     "e.g. bootstraps, bias parameters, etc. will be written.")
    (
     "consistentHits,c",
     po::bool_switch(&(sopt.consistentHits))->default_value(false),
     "Force hits gathered during "
     "quasi-mapping to be \"consistent\" (i.e. co-linear and "
     "approximately the right distance apart).")(
						 "dumpEq", po::bool_switch(&(sopt.dumpEq))->default_value(false),
						 "Dump the equivalence class counts "
						 "that were computed during quasi-mapping")
    ("dumpEqWeights,d",
     po::bool_switch(&(sopt.dumpEqWeights))->default_value(false),
     "Includes \"rich\" equivlance class weights in the output when equivalence "
     "class information is being dumped to file.")
    ("fasterMapping",
     po::bool_switch(&(sopt.fasterMapping))->default_value(false),
     "[Developer]: Disables some extra checks during quasi-mapping. This may make mapping a "
     "little bit faster at the potential cost of returning too many mappings (i.e. some sub-optimal mappings) "
     "for certain reads. Only use this option if you know what it does (enables NIP-skipping)")
    (
     "gcSizeSamp",
     po::value<std::uint32_t>(&(sopt.gcSampFactor))->default_value(1),
     "The value by which to down-sample transcripts when representing the "
     "GC content.  Larger values will reduce memory usage, but may "
     "decrease the fidelity of bias modeling results.")
    (
     "biasSpeedSamp",
     po::value<std::uint32_t>(&(sopt.pdfSampFactor))->default_value(1),
     "The value at which the fragment length PMF is down-sampled "
     "when evaluating sequence-specific & GC fragment bias.  Larger values speed up effective "
     "length correction, but may decrease the fidelity of bias modeling "
     "results.")
    (
     "strictIntersect",
     po::bool_switch(&(sopt.strictIntersect))->default_value(false),
     "Modifies how orphans are "
     "assigned.  When this flag is set, if the intersection of the "
     "quasi-mappings for the left and right "
     "is empty, then all mappings for the left and all mappings for the "
     "right read are reported as orphaned "
     "quasi-mappings")
    (
     "fldMax",
     po::value<size_t>(&(sopt.fragLenDistMax))->default_value(1000),
     "The maximum fragment length to consider when building the empirical "
     "distribution")
    (
     "fldMean",
     po::value<size_t>(&(sopt.fragLenDistPriorMean))->default_value(250),
     "The mean used in the fragment length distribution prior")
    (
     "fldSD",
     po::value<size_t>(&(sopt.fragLenDistPriorSD))->default_value(25),
     "The standard deviation used in the fragment length distribution "
     "prior")
    (
     "forgettingFactor,f",
     po::value<double>(&(sopt.forgettingFactor))->default_value(0.65),
     "The forgetting factor used "
     "in the online learning schedule.  A smaller value results in "
     "quicker learning, but higher variance "
     "and may be unstable.  A larger value results in slower learning but "
     "may be more stable.  Value should "
     "be in the interval (0.5, 1.0].")
    (
     "maxOcc,m", 
     po::value<int>(&(memOptions->max_occ))->default_value(200),
     "(S)MEMs occuring more than this many times won't be considered.")
    (
     "initUniform", po::bool_switch(&(sopt.initUniform))->default_value(false),
     "initialize the offline inference with uniform parameters, rather than seeding with online parameters.")
    (
     "maxReadOcc,w",
     po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(100),
     "Reads \"mapping\" to more than this many places won't be "
     "considered.")
    ("noLengthCorrection",
     po::bool_switch(&(sopt.noLengthCorrection))->default_value(false),
     "[experimental] : Entirely disables length correction when estimating "
     "the abundance of transcripts.  This option can be used with protocols "
     "where one expects that fragments derive from their underlying targets "
     "without regard to that target's length (e.g. QuantSeq)"
     )
    (
     "noEffectiveLengthCorrection",
     po::bool_switch(&(sopt.noEffectiveLengthCorrection))
     ->default_value(false),
     "Disables "
     "effective length correction when computing the "
     "probability that a fragment was generated "
     "from a transcript.  If this flag is passed in, the "
     "fragment length distribution is not taken "
     "into account when computing this probability.")
    (
     "noFragLengthDist",
     po::bool_switch(&(sopt.noFragLengthDist))->default_value(false),
     "[experimental] : "
     "Don't consider concordance with the learned fragment length "
     "distribution when trying to determine "
     "the probability that a fragment has originated from a specified "
     "location.  Normally, Fragments with "
     "unlikely lengths will be assigned a smaller relative probability "
     "than those with more likely "
     "lengths.  When this flag is passed in, the observed fragment length "
     "has no effect on that fragment's "
     "a priori probability.")
    (
     "noBiasLengthThreshold",
     po::bool_switch(&(sopt.noBiasLengthThreshold))->default_value(false),
     "[experimental] : "
     "If this option is enabled, then no (lower) threshold will be set on "
     "how short bias correction can make effective lengths. This can increase the precision "
     "of bias correction, but harm robustness.  The default correction applies a threshold.")
    (
     "numBiasSamples",
     po::value<int32_t>(&numBiasSamples)->default_value(2000000),
     "Number of fragment mappings to use when learning the "
     "sequence-specific bias model.")
    (
     "numAuxModelSamples",
     po::value<uint32_t>(&(sopt.numBurninFrags))->default_value(5000000),
     "The first <numAuxModelSamples> are used to train the "
     "auxiliary model parameters (e.g. fragment length distribution, "
     "bias, etc.).  After ther first <numAuxModelSamples> observations "
     "the auxiliary model parameters will be assumed to have converged "
     "and will be fixed.")
    (
     "numPreAuxModelSamples",
     po::value<uint32_t>(&(sopt.numPreBurninFrags))
     ->default_value(1000000),
     "The first <numPreAuxModelSamples> will have their "
     "assignment likelihoods and contributions to the transcript "
     "abundances computed without applying any auxiliary models.  The "
     "purpose "
     "of ignoring the auxiliary models for the first "
     "<numPreAuxModelSamples> observations is to avoid applying these "
     "models before thier "
     "parameters have been learned sufficiently well.")
    (
     "useVBOpt", po::bool_switch(&(sopt.useVBOpt))->default_value(false),
     "Use the Variational Bayesian EM rather than the "
     "traditional EM algorithm for optimization in the batch passes.")
    (
     "numGibbsSamples",
     po::value<uint32_t>(&(sopt.numGibbsSamples))->default_value(0),
     "Number of Gibbs sampling rounds to "
     "perform.")
    (
     "numBootstraps",
     po::value<uint32_t>(&(sopt.numBootstraps))->default_value(0),
     "Number of bootstrap samples to generate. Note: "
     "This is mutually exclusive with Gibbs sampling.")
    (
     "thinningFactor",
     po::value<uint32_t>(&(sopt.thinningFactor))->default_value(16),
     "Number of steps to discard for every sample kept from the Gibbs chain. "
     "The larger this number, the less chance that subsequent samples are "
     "auto-correlated, but the slower sampling becomes."
    )
    (
     "quiet,q", po::bool_switch(&(sopt.quiet))->default_value(false),
     "Be quiet while doing quantification (don't write informative "
     "output to the console unless something goes wrong).")
    (
     "perTranscriptPrior", po::bool_switch(&(sopt.perTranscriptPrior)), "The "
     "prior (either the default or the argument provided via --vbPrior) will "
     "be interpreted as a transcript-level prior (i.e. each transcript will "
     "be given a prior read count of this value)")
    (
     "vbPrior", po::value<double>(&(sopt.vbPrior))->default_value(1e-3),
     "The prior that will be used in the VBEM algorithm.  This is interpreted "
     "as a per-nucleotide prior, unless the --perTranscriptPrior flag "
     "is also given, in which case this is used as a transcript-level prior")
    (
     "writeOrphanLinks",
     po::bool_switch(&(sopt.writeOrphanLinks))->default_value(false),
     "Write the transcripts that are linked by orphaned reads.")
    (
     "writeUnmappedNames",
     po::bool_switch(&(sopt.writeUnmappedNames))->default_value(false),
     "Write the names of un-mapped reads to the file unmapped_names.txt in the auxiliary directory.")
    ("quasiCoverage,x",
     po::value<double>(&(sopt.quasiCoverage))->default_value(0.0),
     "[Experimental]: The fraction of the read that must be covered by MMPs (of length >= 31) if "
     "this read is to be considered as \"mapped\".  This may help to avoid \"spurious\" mappings. "
     "A value of 0 (the default) denotes no coverage threshold (a single 31-mer can yield a mapping).  "
     "Since coverage by exact matching, large, MMPs is a rather strict condition, this value should likely "
     "be set to something low, if used.");


  po::options_description fmd("\noptions that apply to the old FMD index");
  fmd.add_options()
    (
     "minLen,k",
     po::value<int>(&(memOptions->min_seed_len))->default_value(19),
     "(S)MEMs smaller than this size won't be considered.")
    (
     "sensitive", po::bool_switch(&(sopt.sensitive))->default_value(false),
     "Setting this option enables the splitting of SMEMs that are larger "
     "than 1.5 times the minimum seed length (minLen/k above).  This may "
     "reveal high scoring chains of MEMs "
     "that are masked by long SMEMs.  However, this option makes "
     "lightweight-alignment a bit slower and is "
     "usually not necessary if the reference is of reasonable quality.")
    (
     "extraSensitive",
     po::bool_switch(&(sopt.extraSeedPass))->default_value(false),
     "Setting this option enables an extra pass of \"seed\" search. "
     "Enabling this option may improve sensitivity (the number of reads "
     "having sufficient coverage), but will "
     "typically slow down quantification by ~40%.  Consider enabling this "
     "option if you find the mapping rate to "
     "be significantly lower than expected.")
    (
     "coverage,c", po::value<double>(&coverageThresh)->default_value(0.70),
     "required coverage of read by union of SMEMs to consider it a \"hit\".")
    (
     "splitWidth,s",
     po::value<int>(&(memOptions->split_width))->default_value(0),
     "If (S)MEM occurs fewer than this many times, search for smaller, "
     "contained MEMs. "
     "The default value will not split (S)MEMs, a higher value will "
     "result in more MEMs being explore and, thus, will "
     "result in increased running time.")
    (
     "splitSpanningSeeds,b",
     po::bool_switch(&(sopt.splitSpanningSeeds))->default_value(false),
     "Attempt to split seeds that happen to fall on the "
     "boundary between two transcripts.  This can improve the  fragment "
     "hit-rate, but is usually not necessary.");
    
  po::options_description hidden("\nhidden options");
  hidden.add_options()
    ("numGCBins", po::value<size_t>(&(sopt.numFragGCBins))->default_value(25),
     "Number of bins to use when modeling fragment GC bias")
    (
     "conditionalGCBins", po::value<size_t>(&(sopt.numConditionalGCBins))->default_value(3),
     "Number of different fragment GC models to learn based on read start/end context")
    (
     "numRequiredObs,n",
     po::value(&(sopt.numRequiredFragments))->default_value(50000000),
     "[Deprecated]: The minimum number of observations (mapped reads) "
     "that must be observed before "
     "the inference procedure will terminate.");

  po::options_description testing("\n"
                                  "testing options");
  testing.add_options()
    (
     "posBias", po::value(&(sopt.posBiasCorrect))->zero_tokens(),
     "[experimental] Perform positional bias correction")
    (
      "noRichEqClasses",
      po::bool_switch(&(sopt.noRichEqClasses))->default_value(false),
      "[TESTING OPTION]: Disable \"rich\" equivalent classes.  If this flag is "
      "passed, then "
      "all information about the relative weights for each transcript in the "
      "label of an equivalence class will be ignored, and only the relative "
      "abundance and effective length of each transcript will be considered.")(
      "noFragLenFactor",
      po::bool_switch(&(sopt.noFragLenFactor))->default_value(false),
      "[TESTING OPTION]: Disable the factor in the likelihood that takes into "
      "account the "
      "goodness-of-fit of an alignment with the empirical fragment length "
      "distribution")(
      "rankEqClasses",
      po::bool_switch(&(sopt.rankEqClasses))->default_value(false),
      "[TESTING OPTION]: Keep separate equivalence classes for each distinct "
      "ordering of transcripts in the label.")(
      "noExtrapolateCounts",
      po::bool_switch(&(sopt.dontExtrapolateCounts))->default_value(false),
      "[TESTING OPTION]: When generating posterior counts for Gibbs sampling, "
      "use the directly re-allocated counts in each iteration, rather than extrapolating "
      "from transcript fractions.");

  po::options_description deprecated("\ndeprecated options about which to inform the user");
  deprecated.add_options() (
     "useFSPD", po::bool_switch(&(sopt.useFSPD))->default_value(false),
     "[deprecated] : "
     "Consider / model non-uniformity in the fragment start positions "
     "across the transcript.");

  po::options_description all("salmon quant options");
  all.add(generic).add(advanced).add(testing).add(hidden).add(fmd).add(deprecated);

  po::options_description visible("salmon quant options");
  visible.add(generic).add(advanced);

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
      std::cerr << hstring << std::endl;
      std::cerr << visible << std::endl;
      std::exit(0);
    }

    po::notify(vm);
    
    // If we're supposed to be quiet, set the global logger level to >= warn
    if (sopt.quiet) {
        spdlog::set_level(spdlog::level::warn); //Set global log level to info
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

    ReadExperiment experiment(readLibraries, indexDirectory, sopt);

    // This will be the class in charge of maintaining our
    // rich equivalence classes
    experiment.equivalenceClassBuilder().start();

    auto indexType = experiment.getIndex()->indexType();

    switch (indexType) {
    case SalmonIndexType::FMD: {
      /** Currently no seq-specific bias correction with
       *  FMD index.
       */
      if (sopt.biasCorrect or sopt.gcBiasCorrect) {
        sopt.biasCorrect = false;
        sopt.gcBiasCorrect = false;
        jointLog->warn(
                       "Sequence-specific or fragment GC bias correction require "
                       "use of the quasi-index. Disabling all bias correction");
      }
      quantifyLibrary<SMEMAlignment>(experiment, greedyChain, memOptions, sopt,
                                     coverageThresh, sopt.numThreads);
    } break;
    case SalmonIndexType::QUASI: {
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

      sopt.allowOrphans = true;
      sopt.useQuasi = true;
      quantifyLibrary<QuasiAlignment>(experiment, greedyChain, memOptions, sopt,
                                      coverageThresh, sopt.numThreads);
    } break;
    }

    // Write out information about the command / run
    salmon::utils::writeCmdInfo(sopt, orderedOptions);

    GZipWriter gzw(outputDirectory, jointLog);

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

    free(memOptions);
    size_t tnum{0};

    jointLog->info("writing output \n");

    bfs::path estFilePath = outputDirectory / "quant.sf";

    // Write the main results
    gzw.writeAbundances(sopt, experiment);

    // If we are dumping the equivalence classes, then
    // do it here.
    if (sopt.dumpEq) {
      gzw.writeEquivCounts(sopt, experiment);
    }

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
        //sampler.sampleMultipleChains(experiment, sopt, bsWriter, sopt.numGibbsSamples);
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

    /** If the user requested gene-level abundances, then compute those now **/
    if (vm.count("geneMap")) {
      try {
        salmon::utils::generateGeneLevelEstimates(sopt.geneMapPath,
                                                  outputDirectory);
      } catch (std::invalid_argument& e) {
        fmt::print(stderr, "Error: [{}] when trying to compute gene-level "
                           "estimates. The gene-level file(s) may not exist",
                   e.what());
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
    gzw.writeMeta(sopt, experiment);

  } catch (po::error& e) {
    std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
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
