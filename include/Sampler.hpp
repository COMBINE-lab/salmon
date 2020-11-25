#ifndef __SAMPLER__HPP__
#define __SAMPLER__HPP__

extern "C" {
#include "io_lib/os.h"
#include "io_lib/scram.h"
#undef max
#undef min
}

// for cpp-format
#include "spdlog/fmt/fmt.h"
#include "spdlog/spdlog.h"

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

#include <boost/config.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/program_options.hpp>
#include <boost/timer/timer.hpp>

#include "AlignmentLibrary.hpp"
#include "AlignmentModel.hpp"
#include "BAMQueue.hpp"
#include "ClusterForest.hpp"
#include "FASTAParser.hpp"
#include "FragmentLengthDistribution.hpp"
#include "LibraryFormat.hpp"
#include "MiniBatchInfo.hpp"
#include "OutputUnmappedFilter.hpp"
#include "ReadPair.hpp"
#include "SalmonConfig.hpp"
#include "SalmonMath.hpp"
#include "SalmonOpts.hpp"
#include "SalmonUtils.hpp"
#include "Transcript.hpp"
#include "TranscriptCluster.hpp"

namespace salmon {
namespace sampler {

namespace bfs = boost::filesystem;
using salmon::math::LOG_0;
using salmon::math::LOG_1;
using salmon::math::logAdd;
using salmon::math::logSub;

template <typename FragT> using AlignmentBatch = std::vector<FragT>;

template <typename FragT>
using MiniBatchQueue = tbb::concurrent_queue<MiniBatchInfo<FragT>*>;

template <typename FragT>
using OutputQueue = tbb::concurrent_bounded_queue<FragT*>;

template <typename FragT>
using AlignmentLibraryT = AlignmentLibrary<FragT, EquivalenceClassBuilder<TGValue>>;


template <typename FragT>
void sampleMiniBatch(AlignmentLibraryT<FragT>& alnLib,
                     MiniBatchQueue<AlignmentGroup<FragT*>>& workQueue,
                     std::condition_variable& workAvailable,
                     std::mutex& cvmutex, volatile bool& doneParsing,
                     std::atomic<size_t>& activeBatches,
                     const SalmonOpts& salmonOpts, bool& burnedIn,
                     std::atomic<size_t>& processedReads,
                     OutputQueue<FragT>& outputQueue) {

  // Seed with a real random value, if available
  #if defined(__linux) && defined(__GLIBCXX__) && __GLIBCXX__ >= 20200128
    std::random_device rd("/dev/urandom");
  #else
    std::random_device rd;
  #endif  // defined(__GLIBCXX__) && __GLIBCXX__ >= 2020012
  
  auto log = spdlog::get("jointLog");

  // Create a random uniform distribution
  std::default_random_engine eng(rd());
  std::uniform_real_distribution<> uni(
      0.0, 1.0 + std::numeric_limits<double>::min());

  using salmon::math::LOG_0;
  using salmon::math::LOG_1;
  using salmon::math::LOG_EPSILON;
  using salmon::math::logAdd;
  using salmon::math::logSub;

  bool noLengthCorrection{salmonOpts.noLengthCorrection};
  bool useFragLengthDist{!salmonOpts.noFragLengthDist};
  bool useAuxParams = (processedReads >= salmonOpts.numPreBurninFrags);
  bool considerCondProb = (useAuxParams or burnedIn);

  auto& refs = alnLib.transcripts();
  auto& clusterForest = alnLib.clusterForest();
  auto& fragmentQueue = alnLib.fragmentQueue();
  auto& alignmentGroupQueue = alnLib.alignmentGroupQueue();
  auto& fragLengthDist = *(alnLib.fragmentLengthDistribution());
  auto& alnMod = alnLib.alignmentModel();

  std::vector<FragmentStartPositionDistribution>& fragStartDists =
      alnLib.fragmentStartPositionDistributions();

  const auto expectedLibraryFormat = alnLib.format();

  auto isUnexpectedOrphan = [expectedLibraryFormat](FragT* aln) -> bool {
    return (expectedLibraryFormat.type == ReadType::PAIRED_END and
            !aln->isPaired());
  };

  std::chrono::microseconds sleepTime(1);
  MiniBatchInfo<AlignmentGroup<FragT*>>* miniBatch = nullptr;
  size_t numTranscripts = refs.size();

  while (!doneParsing or !workQueue.empty()) {
    bool foundWork = workQueue.try_pop(miniBatch);
    // If work wasn't immediately available, then wait for it.
    if (!foundWork) {
      std::unique_lock<std::mutex> l(cvmutex);
      workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() {
        return workQueue.try_pop(miniBatch) or doneParsing;
      });
    }
    if (miniBatch != nullptr) {
      ++activeBatches;
      size_t batchReads{0};
      std::vector<AlignmentGroup<FragT*>*>& alignmentGroups =
          *(miniBatch->alignments);

      using TranscriptID = size_t;
      using HitIDVector = std::vector<size_t>;
      using HitProbVector = std::vector<double>;

      std::unordered_map<TranscriptID, std::vector<FragT*>> hitList;
      // Each alignment group corresponds to all of the potential
      // mapping locations of a multi-mapping read
      for (auto alnGroup : alignmentGroups) {
        // Score all of these alignments and sample according to
        // their probabilities
        for (auto a : alnGroup->alignments()) {
          auto transcriptID = a->transcriptID();
          if (transcriptID < 0 or transcriptID >= static_cast<decltype(transcriptID)>(refs.size())) {
            log->warn("Invalid Transcript ID: {}\n", transcriptID);
          }
          hitList[transcriptID].emplace_back(a);
        }
      }

      {
        // Iterate over each group of alignments (a group consists of all
        // alignments reported for a single read).
        for (auto alnGroup : alignmentGroups) {
          double sumOfAlignProbs{LOG_0};
          // update the cluster-level properties
          bool transcriptUnique{true};
          auto firstTranscriptID =
              alnGroup->alignments().front()->transcriptID();
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
            double logFragProb = salmon::math::LOG_1;
            // If we are expecting a paired-end library, and this is an orphan,
            // then logFragProb should be small
            if (isUnexpectedOrphan(aln)) {
              logFragProb = LOG_EPSILON;
            }

            if (flen > 0.0 and aln->isPaired() and useFragLengthDist and
                considerCondProb) {
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

            if (!salmonOpts.noFragLengthDist and aln->fragLen() > 0.0) {
              logFragProb =
                  fragLengthDist.pmf(static_cast<size_t>(aln->fragLen()));
            }

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
            double logAlignCompatProb =
                isCompat ? LOG_1 : salmonOpts.incompatPrior;
            if (!isCompat and salmonOpts.ignoreIncompat) {
              aln->logProb = salmon::math::LOG_0;
              continue;
            }

            // Adjustment to the likelihood due to the
            // error model
            double errLike = salmon::math::LOG_1;
            if (burnedIn and salmonOpts.useErrorModel) {
              errLike = alnMod.logLikelihood(*aln, transcript);
            }

            // Allow for a non-uniform fragment start position distribution
            double startPosProb = -logRefLength;
            auto hitPos = aln->left();

            double auxProb = startPosProb + logFragProb + aln->logQualProb() +
                             errLike + logAlignCompatProb;

            double transcriptLogCount = transcript.mass(false);

            if (transcriptLogCount != LOG_0 and auxProb != LOG_0) {

              aln->logProb = transcriptLogCount + auxProb;
              sumOfAlignProbs = logAdd(sumOfAlignProbs, aln->logProb);

            } else {
              aln->logProb = LOG_0;
            }
          }
          // normalize the hits
          if (sumOfAlignProbs == LOG_0) {
            auto aln = alnGroup->alignments().front();
            log->warn("0 probability fragment [{}] "
                      "encountered\n",
                      aln->getName());
            continue;
          }

          if (transcriptUnique) {
            // avoid r-value ref until we figure out what's
            // up with TBB 4.3
            auto* alnPtr = alnGroup->alignments().front()->clone();
            outputQueue.push(alnPtr);
          } else { // read maps to multiple transcripts
            double r = uni(eng);
            double currentMass{0.0};
            double massInc{0.0};
            bool choseAlignment{false};
            for (auto& aln : alnGroup->alignments()) {
              if (aln->logProb == LOG_0) {
                continue;
              }
              aln->logProb -= sumOfAlignProbs;

              massInc = std::exp(aln->logProb);
              if (currentMass <= r and currentMass + massInc > r) {
                // Write out this read
                // avoid r-value ref until we figure out what's
                // up with TBB 4.3
                auto* alnPtr = aln->clone();
                outputQueue.push(alnPtr);
                currentMass += massInc;
                choseAlignment = true;
                break;
              }
              currentMass += massInc;
            } // end alignment group
            if (BOOST_UNLIKELY(!choseAlignment)) {
              log->warn(
                  "[Sampler.hpp]: Failed to sample an alignment for this read; "
                  "this shouldn't happen\n"
                  "currentMass = {}, r = {}\n",
                  currentMass, r);
            }
          } // non-unique read
          ++batchReads;
        } // end read group
      }   // end timer

      miniBatch->release(fragmentQueue, alignmentGroupQueue);
      delete miniBatch;
      --activeBatches;
      processedReads += batchReads;
    }
    miniBatch = nullptr;
  } // nothing left to process
}

/**
 *  Sample the alignments in the provided library given in current
 *  estimates of transcript abundance.
 */
template <typename FragT>
bool sampleLibrary(AlignmentLibraryT<FragT>& alnLib,
                   const SalmonOpts& salmonOpts, bool burnedIn,
                   bfs::path& sampleFilePath, bool sampleUnaligned) {

  fmt::MemoryWriter msgStr;
  auto log = spdlog::get("jointLog");

  msgStr << "Sampling alignments; outputting results to "
         << sampleFilePath.string() << "\n";

  log->info(msgStr.str());

  auto& refs = alnLib.transcripts();
  size_t numTranscripts = refs.size();
  size_t miniBatchSize{1000};
  size_t numObservedFragments{0};

  MiniBatchQueue<AlignmentGroup<FragT*>> workQueue;
  double logForgettingMass{std::log(1.0)};
  double forgettingFactor{0.60};
  size_t batchNum{0};

  /**
   * Output queue
   */
  volatile bool consumedAllInput{false};
  size_t defaultCapacity = 2000000;
  OutputQueue<FragT> outQueue;
  outQueue.set_capacity(defaultCapacity);

  std::unique_ptr<OutputUnmappedFilter<FragT>> outFilt = nullptr;

  if (sampleUnaligned) {
    outFilt.reset(new OutputUnmappedFilter<FragT>(&outQueue));
  }

  // Reset our reader to the beginning
  if (!alnLib.reset(false, outFilt.get())) {
    fmt::print(stderr,
               "\n\n======== WARNING ========\n"
               "A provided alignment file "
               "is not a regular file and therefore can't be read from "
               "more than once.\n\n"
               "Therefore, we cannot provide sample output alignments "
               "from this file.  No sampled BAM file will be generated. "
               "Please consider re-running Salmon with these alignments "
               "as a regular file!\n"
               "==========================\n\n");

    return false;
  }

  volatile bool doneParsing{false};
  std::condition_variable workAvailable;
  std::mutex cvmutex;
  std::vector<std::thread> workers;
  std::atomic<size_t> activeBatches{0};
  std::atomic<size_t> processedReads{0};

  size_t numProc{0};
  for (uint32_t i = 0; i < salmonOpts.numQuantThreads; ++i) {
    workers.emplace_back(
        sampleMiniBatch<FragT>, std::ref(alnLib), std::ref(workQueue),
        std::ref(workAvailable), std::ref(cvmutex), std::ref(doneParsing),
        std::ref(activeBatches), std::ref(salmonOpts), std::ref(burnedIn),
        std::ref(processedReads), std::ref(outQueue));
  }

  std::thread outputThread(
      [&consumedAllInput, &alnLib, &outQueue, &log, sampleFilePath]() -> void {

        scram_fd* bf = scram_open(sampleFilePath.c_str(), "wb");
        scram_set_option(bf, CRAM_OPT_NTHREADS, 3);
        scram_set_header(bf, alnLib.header());
        scram_write_header(bf);
        if (bf == nullptr) {
          fmt::MemoryWriter errstr;
          errstr << ioutils::SET_RED << "ERROR: " << ioutils::RESET_COLOR
                 << "Couldn't open output bam file " << sampleFilePath.string()
                 << ". Exiting\n";
          log->warn(errstr.str());
          std::exit(-1);
        }

        FragT* aln{nullptr};
        while (!outQueue.empty() or !consumedAllInput) {
          while (outQueue.try_pop(aln)) {
            if (aln != nullptr) {
              int ret = aln->writeToFile(bf);
              if (ret != 0) {
                std::cerr << "ret = " << ret << "\n";
                fmt::MemoryWriter errstr;
                errstr << ioutils::SET_RED << "ERROR:" << ioutils::RESET_COLOR
                       << "Could not write "
                       << "a sampled alignment to the output BAM "
                       << "file. Please check that the file can "
                       << "be created properly and that the disk "
                       << "is not full.  Exiting.\n";
                log->warn(errstr.str());
                std::exit(-1);
              }
              // Eventually, as we do in BAMQueue, we should
              // have queue of bam1_t structures that can be
              // re-used rather than continually calling
              // new and delete.
              delete aln;
              aln = nullptr;
            }
          }
        }

        scram_close(bf); // will delete the header itself
      });

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
      // Don't need to update the batch number or log forgetting mass in this
      // phase
      MiniBatchInfo<AlignmentGroup<FragT*>>* mbi =
          new MiniBatchInfo<AlignmentGroup<FragT*>>(batchNum, alignments,
                                                    logForgettingMass);
      workQueue.push(mbi);
      {
        std::unique_lock<std::mutex> l(cvmutex);
        workAvailable.notify_one();
      }
      alignments = new std::vector<AlignmentGroup<FragT*>*>;
      alignments->reserve(miniBatchSize);
    }
    if (numProc % 1000000 == 0) {
      const char RESET_COLOR[] = "\x1b[0m";
      char green[] = "\x1b[30m";
      green[3] = '0' + static_cast<char>(fmt::GREEN);
      char red[] = "\x1b[30m";
      red[3] = '0' + static_cast<char>(fmt::RED);
      fmt::print(stderr, "\r\r{}processed{} {} {}reads{}", green, red, numProc,
                 green, RESET_COLOR);
    }
    ++numProc;
    alignmentGroupsRemain = bq.getAlignmentGroup(ag);
  }
  std::cerr << "\n";

  // Free the alignments and the vector holding them
  for (auto& aln : *alignments) {
    aln->alignments().clear();
    delete aln;
    aln = nullptr;
  }
  delete alignments;

  doneParsing = true;

  /**
   * This could be a problem for small sets of alignments --- make sure the
   * work queue is empty!!
   * --- Thanks for finding a dataset that exposes this bug, Richard
   * (Smith-Unna)!
   */
  size_t t = 0;
  while (!workQueue.empty()) {
    std::unique_lock<std::mutex> l(cvmutex);
    workAvailable.notify_one();
  }

  size_t tnum{0};
  for (auto& t : workers) {
    fmt::print(stderr, "killing thread {} . . . ", tnum++);
    {
      std::unique_lock<std::mutex> l(cvmutex);
      workAvailable.notify_all();
    }
    t.join();
    fmt::print(stderr, "done\r\r");
  }
  fmt::print(stderr, "\n");
  consumedAllInput = true;

  numObservedFragments += alnLib.numMappedFragments();
  fmt::print(stderr,
             "# observed = {} mapped fragments.\033[F\033[F\033[F\033[F",
             numObservedFragments);

  fmt::print(stderr, "Waiting on output thread\n");
  outputThread.join();
  fmt::print(stderr, "done\n");

  fmt::print(stderr, "\n\n\n\n");
  return true;
}

} // namespace sampler
} // namespace salmon

#endif // __SAMPLER__HPP__
