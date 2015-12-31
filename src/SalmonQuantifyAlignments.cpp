
extern "C" {
#include "io_lib/scram.h"
#include "io_lib/os.h"
}

// for cpp-format
#include "spdlog/details/format.h"

// are these used?
#include <boost/dynamic_bitset.hpp>
#include <boost/lockfree/spsc_queue.hpp>
#include <boost/lockfree/queue.hpp>

#include <tbb/atomic.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <atomic>
#include <vector>
#include <random>
#include <memory>
#include <exception>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <mutex>
#include <thread>
#include <memory>
#include <condition_variable>

#include <tbb/concurrent_queue.h>

#include <boost/timer/timer.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/program_options.hpp>
#include <boost/pending/disjoint_sets.hpp>
#include <boost/math/distributions/normal.hpp>

#include "ClusterForest.hpp"
#include "AlignmentLibrary.hpp"
#include "MiniBatchInfo.hpp"
#include "BAMQueue.hpp"
#include "SalmonMath.hpp"
#include "FASTAParser.hpp"
#include "LibraryFormat.hpp"
#include "Transcript.hpp"
#include "ReadPair.hpp"
#include "ErrorModel.hpp"
#include "AlignmentModel.hpp"
#include "ForgettingMassCalculator.hpp"
#include "FragmentLengthDistribution.hpp"
#include "TranscriptCluster.hpp"
#include "SalmonUtils.hpp"
#include "SalmonConfig.hpp"
#include "SalmonOpts.hpp"
#include "NullFragmentFilter.hpp"
#include "Sampler.hpp"
#include "spdlog/spdlog.h"
#include "EquivalenceClassBuilder.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "GZipWriter.hpp"
#include "TextBootstrapWriter.hpp"

namespace bfs = boost::filesystem;
using salmon::math::LOG_0;
using salmon::math::LOG_1;
using salmon::math::logAdd;
using salmon::math::logSub;

constexpr uint32_t miniBatchSize{1000};

template <typename FragT>
using AlignmentBatch = std::vector<FragT>;

template <typename FragT>
using MiniBatchQueue = tbb::concurrent_queue<MiniBatchInfo<FragT>*>;

using PriorAbundanceVector = std::vector<double>;
using PosteriorAbundanceVector = std::vector<double>;

struct RefSeq {
    RefSeq(char* name, uint32_t len) : RefName(name), RefLength(len) {}
    std::string RefName;
    uint32_t RefLength;
};

/**
 * Multiple each element in the vector `vec` by the factor `scale`.
 */
template <typename T>
void scaleBy(std::vector<T>& vec, T scale) {
    std::for_each(vec.begin(), vec.end(), [scale](T& ele)->void { ele *= scale; });
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
void processMiniBatch(AlignmentLibrary<FragT>& alnLib,
                      ForgettingMassCalculator& fmCalc,
                      uint64_t firstTimestepOfRound,
                      MiniBatchQueue<AlignmentGroup<FragT*>>& workQueue,
                      MiniBatchQueue<AlignmentGroup<FragT*>>* processedCache,
                      std::condition_variable& workAvailable,
                      std::mutex& cvmutex,
                      volatile bool& doneParsing,
                      std::atomic<size_t>& activeBatches,
                      SalmonOpts& salmonOpts,
                      std::atomic<bool>& burnedIn,
                      bool initialRound,
                      std::atomic<size_t>& processedReads) {

    // Seed with a real random value, if available
    std::random_device rd;
    auto& log = salmonOpts.jointLog;

    // Whether or not we are using "banking"
    bool useMassBanking = (!initialRound and salmonOpts.useMassBanking);
    bool useReadCompat = salmonOpts.incompatPrior != salmon::math::LOG_0;

    // Create a random uniform distribution
    std::default_random_engine eng(rd());
    std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());

    //EQClass
    EquivalenceClassBuilder& eqBuilder = alnLib.equivalenceClassBuilder();
    auto& readBias = alnLib.readBias();

    using salmon::math::LOG_0;
    using salmon::math::logAdd;
    using salmon::math::logSub;

    auto& refs = alnLib.transcripts();
    auto& clusterForest = alnLib.clusterForest();
    auto& fragmentQueue = alnLib.fragmentQueue();
    auto& alignmentGroupQueue = alnLib.alignmentGroupQueue();

    std::vector<FragmentStartPositionDistribution>& fragStartDists =
        alnLib.fragmentStartPositionDistributions();

    auto& fragLengthDist = *(alnLib.fragmentLengthDistribution());
    auto& alnMod = alnLib.alignmentModel();

    bool useFSPD{!salmonOpts.noFragStartPosDist};
    bool useFragLengthDist{!salmonOpts.noFragLengthDist};
    bool noFragLenFactor{salmonOpts.noFragLenFactor};

    double startingCumulativeMass = fmCalc.cumulativeLogMassAt(firstTimestepOfRound);
    const auto expectedLibraryFormat = alnLib.format();
    uint32_t numBurninFrags{salmonOpts.numBurninFrags};

    bool useAuxParams = (processedReads > salmonOpts.numPreBurninFrags);

    std::chrono::microseconds sleepTime(1);
    MiniBatchInfo<AlignmentGroup<FragT*>>* miniBatch = nullptr;
    bool updateCounts = initialRound;
    size_t numTranscripts = refs.size();

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
            workAvailable.wait(l, [&miniBatch, &workQueue, &doneParsing]() { return workQueue.try_pop(miniBatch) or doneParsing; });
        }

	    // If we actually got some work
        if (miniBatch != nullptr) {

            useAuxParams = (processedReads > salmonOpts.numPreBurninFrags);
            ++activeBatches;
            size_t batchReads{0};

            // double logForgettingMass = fmCalc();
            double logForgettingMass{0.0};
            uint64_t currentMinibatchTimestep{0};
	    // logForgettingMass and currentMinibatchTimestep are OUT parameters!
            fmCalc.getLogMassAndTimestep(logForgettingMass, currentMinibatchTimestep);
            miniBatch->logForgettingMass = logForgettingMass;

            std::vector<AlignmentGroup<FragT*>*>& alignmentGroups = *(miniBatch->alignments);

            using TranscriptID = size_t;
            using HitIDVector = std::vector<size_t>;
            using HitProbVector = std::vector<double>;

	    {
                // Iterate over each group of alignments (a group consists of all alignments reported
                // for a single read).  Distribute the read's mass proportionally dependent on the
                // current
                for (auto alnGroup : alignmentGroups) {

                    // EQCLASS
                    std::vector<uint32_t> txpIDs;
                    std::vector<double> auxProbs;
                    double auxDenom = salmon::math::LOG_0;

		    // The alignments must be sorted by transcript id
                    alnGroup->sortHits();

                    double sumOfAlignProbs{LOG_0};

                    // update the cluster-level properties
                    bool transcriptUnique{true};
                    auto firstTranscriptID = alnGroup->alignments().front()->transcriptID();
                    std::unordered_set<size_t> observedTranscripts;

                    for (auto& aln : alnGroup->alignments()) {
                        auto transcriptID = aln->transcriptID();
                        auto& transcript = refs[transcriptID];
                        transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                        double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;
                        double logFragProb = salmon::math::LOG_1;

                        if (!salmonOpts.noFragLengthDist and useAuxParams) {
                            /** Forget reads that are not paired **/
                            /*
                            if(aln->fragLen() == 0) {
                                if (aln->isLeft() and transcript.RefLength - aln->left() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(transcript.RefLength - aln->left());
                                } else if (aln->isRight() and aln->right() < fragLengthDist.maxVal()) {
                                    logFragProb = fragLengthDist.cmf(aln->right());
                                }
                            } else {
                            }
                            */
                            if(aln->isPaired() and aln->fragLen() > 0) {
                                logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln->fragLen()));
                            }
                        }

                        // TESTING
                        if (noFragLenFactor) { logFragProb = LOG_1; }

                        // @TODO: handle this case better
                        //double fragProb = cdf(fragLengthDist, fragLength + 0.5) - cdf(fragLengthDist, fragLength - 0.5);
                        //fragProb = std::max(fragProb, 1e-3);
                        //fragProb /= cdf(fragLengthDist, refLength);

                        // The alignment probability is the product of a transcript-level term (based on abundance and) an alignment-level
                        // term below which is P(Q_1) * P(Q_2) * P(F | T)
                        double logRefLength = std::log(refLength);

                        // The probability that the fragments align to the given strands in the
                        // given orientations.
                        double logAlignCompatProb = (useReadCompat) ?
                                (salmon::utils::logAlignFormatProb(aln->libFormat(), expectedLibraryFormat, salmonOpts.incompatPrior)) :
                                LOG_1;

                        // Adjustment to the likelihood due to the
                        // error model
                        double errLike = salmon::math::LOG_1;
                        //if (burnedIn and salmonOpts.useErrorModel) {
                        if (useAuxParams and salmonOpts.useErrorModel) {
                            errLike = alnMod.logLikelihood(*aln, transcript);
                        }

                        // Allow for a non-uniform fragment start position distribution
                        double startPosProb = -logRefLength;
                        auto hitPos = aln->left();
                        if (useFSPD and burnedIn and hitPos < refLength) {
                            auto& fragStartDist =
                                fragStartDists[transcript.lengthClassIndex()];
                            startPosProb = fragStartDist(hitPos, refLength, logRefLength);
                        }

                        // The total auxiliary probabilty is the product (sum in log-space) of
                        // The fragment length probabilty
                        // The mapping score (under error model) probability
                        // The fragment compatibility probability

                        // The auxProb does *not* account for the start position
                        // probability!
                        double auxProb = logFragProb + errLike + logAlignCompatProb;

                        // The overall mass of this transcript, which is used to
                        // account for this transcript's relaive abundance
                        double transcriptLogCount = transcript.mass(initialRound);

                        if ( transcriptLogCount != LOG_0 and
                              auxProb != LOG_0 and
                              startPosProb != LOG_0 ) {
                            aln->logProb = transcriptLogCount + auxProb + startPosProb;

                            sumOfAlignProbs = logAdd(sumOfAlignProbs, aln->logProb);
                            if (updateCounts and
                                    observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
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
                        auto aln = alnGroup->alignments().front();
                        log->warn("0 probability fragment [{}] "
                                  "encountered \n", aln->getName());
                        continue;
                    }

                    // EQCLASS
                    double auxProbSum{0.0};
                    for (auto& p : auxProbs) {
                        p = std::exp(p - auxDenom);
                        auxProbSum += p;
                    }

                    if (txpIDs.size() > 0) {
                        TranscriptGroup tg(txpIDs);
                        eqBuilder.addGroup(std::move(tg), auxProbs);
                    }


                    // Are we doing bias correction?
                    bool needBiasSample = salmonOpts.biasCorrect;

                    // Normalize the scores
                    for (auto& aln : alnGroup->alignments()) {
                        if (aln->logProb == LOG_0) { continue; }
                        aln->logProb -= sumOfAlignProbs;

                        auto transcriptID = aln->transcriptID();
                        auto& transcript = refs[transcriptID];

                        double newMass = logForgettingMass + aln->logProb;
                        transcript.addMass(newMass);
                        transcript.setLastTimestepUpdated(currentMinibatchTimestep);

                        /**
                         * Update the auxiliary models.
                         **/
                        double r = uni(eng);
                        if (!burnedIn and r < std::exp(aln->logProb)) {
                            /**
                             * Update the bias sequence-specific bias model
                             **/
                            if (needBiasSample and salmonOpts.numBiasSamples > 0) {
                                // the "start" position is the leftmost position if
                                // we hit the forward strand, and the leftmost
                                // position + the read length if we hit the reverse complement
                                bam_seq_t* r = aln->get5PrimeRead();
                                if (r) {
                                    bool fwd{bam_strand(r) == 0};
                                    int32_t pos{bam_pos(r)};
                                    int32_t startPos = fwd ? pos : pos + bam_seq_len(r);
                                    auto dir = salmon::utils::boolToDirection(fwd);

                                    if (startPos > 0 and startPos < transcript.RefLength) {
                                        const char* txpStart = transcript.Sequence;
                                        const char* readStart = txpStart + startPos;
                                        const char* txpEnd = txpStart + transcript.RefLength;
                                        bool success = readBias.update(txpStart, readStart, txpEnd, dir);
                                        if (success) {
                                            salmonOpts.numBiasSamples -= 1;
                                            needBiasSample = false;
                                        }
                                    }
                                }
                            }

                            // Update the error model
                            if (salmonOpts.useErrorModel) {
                                alnMod.update(*aln, transcript, LOG_1, logForgettingMass);
                            }
                            // Update the fragment length distribution
                            if (aln->isPaired() and !salmonOpts.noFragLengthDist) {
                                double fragLength = aln->fragLen();
                                fragLengthDist.addVal(fragLength, logForgettingMass);
                            }
                            // Update the fragment start position distribution
                            if (useFSPD) {
                                auto hitPos = aln->left();
                                auto& fragStartDist =
                                    fragStartDists[transcript.lengthClassIndex()];
                                fragStartDist.addVal(hitPos,
                                        transcript.RefLength,
                                        logForgettingMass);
                            }
                        }
                    }

                    // update the single target transcript
                    if (transcriptUnique) {
                        if (updateCounts) {
                            refs[firstTranscriptID].addUniqueCount(1);
                        }
                        clusterForest.updateCluster(firstTranscriptID, 1,
                                                    logForgettingMass, updateCounts);
                    } else { // or the appropriate clusters
                        // ughh . . . C++ still has some very rough edges
                        clusterForest.template mergeClusters<FragT>(alnGroup->alignments().begin(),
                                                           alnGroup->alignments().end());
                        clusterForest.updateCluster(alnGroup->alignments().front()->transcriptID(),
                                                    1, logForgettingMass, updateCounts);
                    }

                    ++batchReads;
                } // end read group
            }// end timer

            double individualTotal = LOG_0;
            {
                /*
                // M-step
                for (auto kv = hitList.begin(); kv != hitList.end(); ++kv) {
                    auto transcriptID = kv->first;
                    // The target must be a valid transcript
                    if (transcriptID >= numTranscripts or transcriptID < 0) {std::cerr << "index " << transcriptID << " out of bounds\n"; }

                    auto& transcript = refs[transcriptID];

                    // The prior probability
                    double hitMass{LOG_0};

                    // The set of alignments that match transcriptID
                    auto& hits = kv->second;
                    std::for_each(hits.begin(), hits.end(), [&](FragT* aln) -> void {
                            if (!std::isfinite(aln->logProb)) { log->warn("hitMass = {}\n", aln->logProb); }
                            hitMass = logAdd(hitMass, aln->logProb);
                    });

                    // Lock the target
                    if (hitMass == LOG_0) {
                        log->warn("\n\n\n\nA set of *valid* alignments for a read appeared to "
                                  "have 0 probability.  This should not happen.  Please report "
                                  "this bug.  exiting!");

                        std::cerr << "\n\n\n\nA set of *valid* alignments for a read appeared to "
                                  << "have 0 probability.  This should not happen.  Please report "
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
                if (useFSPD) {
                    // update all of the fragment start position
                    // distributions
                    for (auto& fspd : fragStartDists) {
                        fspd.update();
                    }
                }
                fragLengthDist.cacheCMF();
                // NOTE: only one thread should succeed here, and that
                // thread will set burnedIn to true
                alnLib.updateTranscriptLengthsAtomic(burnedIn);
            }
        }
        miniBatch = nullptr;
    } // nothing left to process
}

/**
  *  Quantify the targets given in the file `transcriptFile` using the
  *  alignments given in the file `alignmentFile`, and write the results
  *  to the file `outputFile`.  The reads are assumed to be in the format
  *  specified by `libFmt`.
  *
  */
template <typename FragT>
bool quantifyLibrary(
        AlignmentLibrary<FragT>& alnLib,
        size_t numRequiredFragments,
        SalmonOpts& salmonOpts) {

    std::atomic<bool> burnedIn{false};

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
    size_t maxCacheSize{salmonOpts.mappingCacheMemoryLimit};

    NullFragmentFilter<FragT>* nff = nullptr;
    bool terminate{false};

    // Give ourselves some space
    fmt::print(stderr, "\n\n\n\n");

    while (numObservedFragments < numRequiredFragments and !terminate) {
        if (!initialRound) {

    	    size_t numToCache = (useMassBanking) ?
				(alnLib.numMappedFragments() - alnLib.numUniquelyMappedFragments()) :
				(alnLib.numMappedFragments());

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
	        		true, 		/* increment # of passes */
		        	nff, 		/* fragment filter */
        			useMassBanking  /* only process ambiguously mapped fragments*/
	        	)) {
                fmt::print(stderr,
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
        auto currentQuantThreads = (haveCache) ?
                                   salmonOpts.numQuantThreads + salmonOpts.numParseThreads :
                                   salmonOpts.numQuantThreads;

        uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();
        if (firstTimestepOfRound > 1) {
            firstTimestepOfRound -= 1;
        }

        for (uint32_t i = 0; i < currentQuantThreads; ++i) {
            workers.emplace_back(processMiniBatch<FragT>,
                    std::ref(alnLib),
                    std::ref(fmCalc),
                    firstTimestepOfRound,
                    std::ref(*workQueuePtr),
                    processedCachePtr,
                    std::ref(workAvailable), std::ref(cvmutex),
                    std::ref(doneParsing), std::ref(activeBatches),
                    std::ref(salmonOpts),
                    std::ref(burnedIn),
                    initialRound,
                    std::ref(totalProcessedReads));
        }

        if (!haveCache) {
            size_t numProc{0};

            BAMQueue<FragT>& bq = alnLib.getAlignmentGroupQueue();
            std::vector<AlignmentGroup<FragT*>*>* alignments = new std::vector<AlignmentGroup<FragT*>*>;
            alignments->reserve(miniBatchSize);
            AlignmentGroup<FragT*>* ag;

            bool alignmentGroupsRemain = bq.getAlignmentGroup(ag);
            while (alignmentGroupsRemain or alignments->size() > 0) {
                if (alignmentGroupsRemain) { alignments->push_back(ag); }
                // If this minibatch has reached the size limit, or we have nothing
                // left to fill it up with
                if (alignments->size() >= miniBatchSize or !alignmentGroupsRemain) {
                    ++batchNum;
                    double logForgettingMass = 0.0;
                    MiniBatchInfo<AlignmentGroup<FragT*>>* mbi =
                        new MiniBatchInfo<AlignmentGroup<FragT*>>(batchNum, alignments, logForgettingMass);
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
                    delete alnGroup; alnGroup = nullptr;
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
          * --- Thanks for finding a dataset that exposes this bug, Richard (Smith-Unna)!
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

        fmt::print(stderr, "# observed = {} / # required = {}\033[A\033[A\033[A\033[A\033[A",
                   numObservedFragments, numRequiredFragments);

        if (initialRound) {
            salmonOpts.jointLog->info("\n\n\nCompleted first pass through the alignment file.\n"
                                      "Total # of mapped reads : {}\n"
                                      "# of uniquely mapped reads : {}\n"
                                      "# ambiguously mapped reads : {}\n\n\n",
                                      alnLib.numMappedFragments(),
                                      alnLib.numUniquelyMappedFragments(),
                                      alnLib.numMappedFragments() -
                                      alnLib.numUniquelyMappedFragments());
        }

        initialRound = false;

        // If we're done our second pass and we've decided to use
        // the in-memory cache, then activate it now.
        if (!initialRound and processedCachePtr != nullptr) {
            haveCache = true;
        }
        //EQCLASS
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
        salmonOpts.jointLog->warn("Only {} fragments were mapped, but the number of burn-in fragments was set to {}.\n"
                "The effective lengths have been computed using the observed mappings.\n",
                alnLib.numMappedFragments(), salmonOpts.numBurninFrags);
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
bool processSample(AlignmentLibrary<ReadT>& alnLib,
                   const std::string& runStartTime,
                   size_t requiredObservations,
                   SalmonOpts& sopt,
                   boost::filesystem::path outputDirectory) {

    auto& jointLog = sopt.jointLog;
    // EQCLASS
    alnLib.equivalenceClassBuilder().start();

    bool burnedIn = quantifyLibrary<ReadT>(alnLib, requiredObservations, sopt);

    // EQCLASS
    // NOTE: A side-effect of calling the optimizer is that
    // the `EffectiveLength` field of each transcript is
    // set to its final value.
    CollapsedEMOptimizer optimizer;
    jointLog->info("starting optimizer");
    salmon::utils::normalizeAlphas(sopt, alnLib);
    bool optSuccess = optimizer.optimize(alnLib, sopt, 0.01, 10000);
    // If the optimizer didn't work, then bail out here.
    if (!optSuccess) { return false; }
    jointLog->info("finished optimizer");

    // EQCLASS
    fmt::print(stderr, "\n\nwriting output \n");
    GZipWriter gzw(outputDirectory, jointLog);
    // Write the main results
    gzw.writeAbundances(sopt, alnLib);
    // Write meta-information about the run
    gzw.writeMeta(sopt, alnLib, runStartTime);

    if (sopt.numGibbsSamples > 0) {

        jointLog->info("Starting Gibbs Sampler");
        CollapsedGibbsSampler sampler;
        // The function we'll use as a callback to write samples
        std::function<bool(const std::vector<int>&)> bsWriter =
            [&gzw](const std::vector<int>& alphas) -> bool {
                return gzw.writeBootstrap(alphas);
            };

        bool sampleSuccess = sampler.sample(alnLib, sopt,
                bsWriter,
                sopt.numGibbsSamples);
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
        bool bootstrapSuccess = optimizer.gatherBootstraps(
                alnLib, sopt,
                bsWriter, 0.01, 10000);
        jointLog->info("Finished Bootstrapping");
        if (!bootstrapSuccess) {
            jointLog->error("Encountered error during bootstrapping.\n"
                    "This should not happen.\n"
                    "Please file a bug report on GitHub.\n");
            return false;
        }
    }



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
        bool didSample = salmon::sampler::sampleLibrary<ReadT>(alnLib, sopt, burnedIn, sampleFilePath, sopt.sampleUnaligned);
        if (!didSample) {
            jointLog->warn("There may have been a problem generating the sampled output file; please check the log\n");
        }
    }

    return true;
}


int salmonAlignmentQuantify(int argc, char* argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace po = boost::program_options;
    namespace bfs = boost::filesystem;

    SalmonOpts sopt;

    uint32_t numThreads{4};
    size_t requiredObservations{50000000};

    po::options_description basic("\nbasic options");
    basic.add_options()
    ("version,v", "print version string.")
    ("help,h", "produce help message.")
    ("libType,l", po::value<std::string>()->required(), "Format string describing the library type.")
    ("alignments,a", po::value<vector<string>>()->multitoken()->required(), "input alignment (BAM) file(s).")
    ("targets,t", po::value<std::string>()->required(), "FASTA format file containing target transcripts.")
    ("threads,p", po::value<uint32_t>(&numThreads)->default_value(6), "The number of threads to use concurrently. "
                                            "The alignment-based quantification mode of salmon is usually I/O bound "
                                            "so until there is a faster multi-threaded SAM/BAM parser to feed the "
                                            "quantification threads, one should not expect much of a speed-up beyond "
                                            "~6 threads.")
    ("biasCorrect", po::value(&(sopt.biasCorrect))->zero_tokens(), "Perform sequence-specific bias correction.")
    ("incompatPrior", po::value<double>(&(sopt.incompatPrior))->default_value(1e-20), "This option "
                        "sets the prior probability that an alignment that disagrees with the specified "
                        "library type (--libType) results from the true fragment origin.  Setting this to 0 "
                        "specifies that alignments that disagree with the library type should be \"impossible\", "
                        "while setting it to 1 says that alignments that disagree with the library type are no "
                        "less likely than those that do")
    ("useErrorModel", po::bool_switch(&(sopt.useErrorModel))->default_value(false), "[Currently Experimental] : "
                        "Learn and apply an error model for the aligned reads.  This takes into account the "
                        "the observed frequency of different types of mismatches when computing the likelihood of "
                        "a given alignment.")
    ("output,o", po::value<std::string>()->required(), "Output quantification directory.")
    ("numRequiredObs,n", po::value(&requiredObservations)->default_value(50000000),
                                        "[Deprecated]: The minimum number of observations (mapped reads) that must be observed before "
                                        "the inference procedure will terminate.  If fewer mapped reads exist in the "
                                        "input file, then it will be read through multiple times.")
    ("geneMap,g", po::value<std::string>(), "File containing a mapping of transcripts to genes.  If this file is provided "
                                        "Salmon will output both quant.sf and quant.genes.sf files, where the latter "
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping "
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format "
                                        "where each line contains the name of a transcript and the gene to which it belongs "
                                        "separated by a tab.  The extension of the file is used to determine how the file "
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF "
                                        "format; files with any other extension are assumed to be in the simple format.");

    // no sequence bias for now
    sopt.useMassBanking = false;
    sopt.noSeqBiasModel = true;
    sopt.noRichEqClasses = false;

    po::options_description advanced("\nadvanced options");
    advanced.add_options()
    ("fldMax" , po::value<size_t>(&(sopt.fragLenDistMax))->default_value(800), "The maximum fragment length to consider when building the empirical distribution")
    ("fldMean", po::value<size_t>(&(sopt.fragLenDistPriorMean))->default_value(200), "The mean used in the fragment length distribution prior")
    ("fldSD" , po::value<size_t>(&(sopt.fragLenDistPriorSD))->default_value(80), "The standard deviation used in the fragment length distribution prior")
    ("forgettingFactor,f", po::value<double>(&(sopt.forgettingFactor))->default_value(0.65), "The forgetting factor used "
                        "in the online learning schedule.  A smaller value results in quicker learning, but higher variance "
                        "and may be unstable.  A larger value results in slower learning but may be more stable.  Value should "
                        "be in the interval (0.5, 1.0].")
    ("mappingCacheMemoryLimit", po::value<uint32_t>(&(sopt.mappingCacheMemoryLimit))->default_value(2000000), "If the file contained fewer than this "
                                        "many mapped reads, then just keep the data in memory for subsequent rounds of inference. Obviously, this value should "
                                        "not be too large if you wish to keep a low memory usage, but setting it large enough to accommodate all of the mapped "
                                        "read can substantially speed up inference on \"small\" files that contain only a few million reads.")
    ("maxReadOcc,w", po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(200), "Reads \"mapping\" to more than this many places won't be considered.")
    ("noEffectiveLengthCorrection", po::bool_switch(&(sopt.noEffectiveLengthCorrection))->default_value(false), "Disables "
                        "effective length correction when computing the probability that a fragment was generated "
                        "from a transcript.  If this flag is passed in, the fragment length distribution is not taken "
                        "into account when computing this probability.")
    ("noFragLengthDist", po::bool_switch(&(sopt.noFragLengthDist))->default_value(false), "[Currently Experimental] : "
                        "Don't consider concordance with the learned fragment length distribution when trying to determine "
                        "the probability that a fragment has originated from a specified location.  Normally, Fragments with "
                         "unlikely lengths will be assigned a smaller relative probability than those with more likely "
                        "lengths.  When this flag is passed in, the observed fragment length has no effect on that fragment's "
                        "a priori probability.")
    ("noFragStartPosDist", po::bool_switch(&(sopt.noFragStartPosDist))->default_value(false), "[Currently Experimental] : "
                        "Don't consider / model non-uniformity in the fragment start positions "
                        "across the transcript.")
    /*
    // Don't expose this yet
    ("noRichEqClasses", po::bool_switch(&(sopt.noRichEqClasses))->default_value(false),
                        "Disable \"rich\" equivalent classes.  If this flag is passed, then "
                        "all information about the relative weights for each transcript in the "
                        "label of an equivalence class will be ignored, and only the relative "
                        "abundance and effective length of each transcript will be considered.")
                        */
    ("numErrorBins", po::value<uint32_t>(&(sopt.numErrorBins))->default_value(6), "The number of bins into which to divide "
                        "each read when learning and applying the error model.  For example, a value of 10 would mean that "
                        "effectively, a separate error model is leared and applied to each 10th of the read, while a value of "
                        "3 would mean that a separate error model is applied to the read beginning (first third), middle (second third) "
                        "and end (final third).")
    ("numPreAuxModelSamples", po::value<uint32_t>(&(sopt.numPreBurninFrags))->default_value(1000000), "The first <numPreAuxModelSamples> will have their "
     			"assignment likelihoods and contributions to the transcript abundances computed without applying any auxiliary models.  The purpose "
			"of ignoring the auxiliary models for the first <numPreAuxModelSamples> observations is to avoid applying these models before thier "
			"parameters have been learned sufficiently well.")
    ("numAuxModelSamples", po::value<uint32_t>(&(sopt.numBurninFrags))->default_value(5000000), "The first <numAuxModelSamples> are used to train the "
     			"auxiliary model parameters (e.g. fragment length distribution, bias, etc.).  After ther first <numAuxModelSamples> observations "
			"the auxiliary model parameters will be assumed to have converged and will be fixed.")
    ("sampleOut,s", po::bool_switch(&(sopt.sampleOutput))->default_value(false), "Write a \"postSample.bam\" file in the output directory "
                        "that will sample the input alignments according to the estimated transcript abundances. If you're "
                        "going to perform downstream analysis of the alignments with tools which don't, themselves, take "
                        "fragment assignment ambiguity into account, you should use this output.")
    ("sampleUnaligned,u", po::bool_switch(&(sopt.sampleUnaligned))->default_value(false), "In addition to sampling the aligned reads, also write "
                        "the un-aligned reads to \"posSample.bam\".")
    ("numGibbsSamples", po::value<uint32_t>(&(sopt.numGibbsSamples))->default_value(0), "[*super*-experimental]: Number of Gibbs sampling rounds to "
     "perform.")
    ("numBootstraps", po::value<uint32_t>(&(sopt.numBootstraps))->default_value(0), "[experimental]: Number of bootstrap samples to generate. Note: "
      "This is mutually exclusive with Gibbs sampling.");

    po::options_description testing("\n"
            "testing options");
    testing.add_options()
        ("noRichEqClasses", po::bool_switch(&(sopt.noRichEqClasses))->default_value(false),
                        "[TESTING OPTION]: Disable \"rich\" equivalent classes.  If this flag is passed, then "
                        "all information about the relative weights for each transcript in the "
                        "label of an equivalence class will be ignored, and only the relative "
                        "abundance and effective length of each transcript will be considered.")
        ("noFragLenFactor", po::bool_switch(&(sopt.noFragLenFactor))->default_value(false),
                        "[TESTING OPTION]: Disable the factor in the likelihood that takes into account the "
                        "goodness-of-fit of an alignment with the empirical fragment length "
                        "distribution");


    po::options_description all("salmon quant options");
    all.add(basic).add(advanced).add(testing);

    po::options_description visible("salmon quant options");
    visible.add(basic).add(advanced);

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(all).run();

        po::store(orderedOptions, vm);

        if (vm.count("help")) {
            std::cout << "Salmon quant (alignment-based)\n";
            std::cout << visible << std::endl;
            std::exit(0);
        }
        po::notify(vm);

        sopt.alnMode = true;

        if (numThreads < 2) {
            fmt::print(stderr, "salmon requires at least 2 threads --- "
                               "setting # of threads = 2\n");
            numThreads = 2;
        }
        sopt.numThreads = numThreads;

        if (sopt.forgettingFactor <= 0.5 or
            sopt.forgettingFactor > 1.0) {
            fmt::print(stderr, "The forgetting factor must be in (0.5, 1.0], "
                               "but the value {} was provided\n", sopt.forgettingFactor);
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
        fmt::print(stderr, "{}", commentString);

        // TODO: Fix fragment start pos dist
        // sopt.noFragStartPosDist = true;

        // Get the time at the start of the run
        std::time_t result = std::time(NULL);
        std::string runStartTime(std::asctime(std::localtime(&result)));
        runStartTime.pop_back(); // remove the newline

        // Verify the geneMap before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("geneMap")) {
            // Make sure the provided file exists
            geneMapPath = vm["geneMap"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                fmt::print(stderr, "Could not find transcript <=> gene map file {} \n"
                           "Exiting now; please either omit the \'geneMap\' option or "
                           "provide a valid file\n", geneMapPath);
                std::exit(1);
            }
        }

        vector<string> alignmentFileNames = vm["alignments"].as<vector<string>>();
        vector<bfs::path> alignmentFiles;
        for (auto& alignmentFileName : alignmentFileNames) {
            bfs::path alignmentFile(alignmentFileName);
            if (!bfs::exists(alignmentFile)) {
                std::stringstream ss;
                ss << "The provided alignment file: " << alignmentFile <<
                    " does not exist!\n";
                throw std::invalid_argument(ss.str());
            } else {
                alignmentFiles.push_back(alignmentFile);
            }
        }

        std::string libFmtStr = vm["libType"].as<std::string>();
        LibraryFormat libFmt = salmon::utils::parseLibraryFormatStringNew(libFmtStr);
        if (libFmt.check()) {
            std::cerr << libFmt << "\n";
        } else {
            std::stringstream ss;
            ss << libFmt << " is invalid!";
            throw std::invalid_argument(ss.str());
        }

        bfs::path outputDirectory(vm["output"].as<std::string>());
        // If the path exists
        if (bfs::exists(outputDirectory)) {
            // If it is not a directory, then complain
            if (!bfs::is_directory(outputDirectory)) {
                std::stringstream errstr;
                errstr << "Path [" << outputDirectory << "] already exists "
                       << "and is not a directory.\n"
                       << "Please either remove this file or choose another "
                       << "output path.\n";
                throw std::invalid_argument(errstr.str());
            }
        } else { // If the path doesn't exist, then create it
            if (!bfs::create_directories(outputDirectory)) { // creation failed for some reason
                std::stringstream errstr;
                errstr << "Could not create output directory ["
                       << outputDirectory << "], please check that it is valid.";
                throw std::invalid_argument(errstr.str());
            }
        }

        bfs::path logDirectory = outputDirectory / "logs";

        // Create the logger and the logging directory
        bfs::create_directories(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "Logs will be written to " << logDirectory.string() << "\n";

        bfs::path logPath = logDirectory / "salmon.log";
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("consoleLog", {consoleSink});
        auto fileLog = spdlog::create("fileLog", {fileSink});
        auto jointLog = spdlog::create("jointLog", {fileSink, consoleSink});

        sopt.jointLog = jointLog;
        sopt.fileLog = fileLog;

        // Verify that no inconsistent options were provided
        if (sopt.numGibbsSamples > 0 and sopt.numBootstraps > 0) {
            jointLog->error("You cannot perform both Gibbs sampling and bootstrapping. "
                            "Please choose one.");
            jointLog->flush();
            std::exit(1);
        }

        if (!sopt.sampleOutput and sopt.sampleUnaligned) {
            fmt::MemoryWriter wstr;
            wstr << "WARNING: you passed in the (-u/--sampleUnaligned) flag, but did not request a sampled "
                 << "output file (-s/--sampleOut).  This flag will be ignored!\n";
            jointLog->warn() << wstr.str();
        }

        // maybe arbitrary, but if it's smaller than this, consider it
        // equal to LOG_0
        if (sopt.incompatPrior < 1e-320) {
            sopt.incompatPrior = salmon::math::LOG_0;
        } else {
            sopt.incompatPrior = std::log(sopt.incompatPrior);
        }

        // Now create a subdirectory for any parameters of interest
        bfs::path paramsDir = outputDirectory / "libParams";
        if (!boost::filesystem::exists(paramsDir)) {
            if (!boost::filesystem::create_directories(paramsDir)) {
                fmt::print(stderr, "{}ERROR{}: Could not create "
                           "output directory for experimental parameter "
                           "estimates [{}]. exiting.", ioutils::SET_RED,
                           ioutils::RESET_COLOR, paramsDir);
                std::exit(-1);
            }
        }


        // Write out information about the command / run
        {
            bfs::path cmdInfoPath = outputDirectory / "cmd_info.json";
            std::ofstream os(cmdInfoPath.string());
            cereal::JSONOutputArchive oa(os);
            oa(cereal::make_nvp("salmon_version", std::string(salmon::version)));
            for (auto& opt : orderedOptions.options) {
                if (opt.value.size() == 1) {
                    oa(cereal::make_nvp(opt.string_key, opt.value.front()));
                } else {
                    oa(cereal::make_nvp(opt.string_key, opt.value));
                }
            }
        }

        // The transcript file contains the target sequences
        bfs::path transcriptFile(vm["targets"].as<std::string>());

        // Currently, one thread is used for parsing the alignment file.
        // Hopefully, in the future, samtools will implemented multi-threaded
        // BAM/SAM parsing, as this is the current bottleneck.  For the time
        // being, however, the number of quantification threads is the
        // total number of threads - 1.
        uint32_t numParseThreads = std::min(uint32_t(6),
                                            std::max(uint32_t(2), uint32_t(std::ceil(numThreads/2.0))));
        numThreads = std::max(numThreads, numParseThreads);
        uint32_t numQuantThreads = std::max(uint32_t(2), uint32_t(numThreads - numParseThreads));
        sopt.numQuantThreads = numQuantThreads;
        sopt.numParseThreads = numParseThreads;
        std::cerr << "numQuantThreads = " << numQuantThreads << "\n";

        bool success{false};

        switch (libFmt.type) {
            case ReadType::SINGLE_END:
                {
                    AlignmentLibrary<UnpairedRead> alnLib(alignmentFiles,
                                                          transcriptFile,
                                                          libFmt,
                                                          sopt);

                    success = processSample<UnpairedRead>(alnLib, runStartTime,
                                                          requiredObservations, sopt,
                                                          outputDirectory);
                }
                break;
            case ReadType::PAIRED_END:
                {
                    AlignmentLibrary<ReadPair> alnLib(alignmentFiles,
                                                      transcriptFile,
                                                      libFmt,
                                                      sopt);

                    success = processSample<ReadPair>(alnLib, runStartTime,
                                                      requiredObservations, sopt,
                                                      outputDirectory);
                }
                break;
            default:
                std::cerr << "Cannot quantify library of unknown format "
                          << libFmt << "\n";
                std::exit(1);
        }

        // Make sure the quantification was successful.
        if (!success) {
            jointLog->error("Quantification was un-successful.  Please check the log "
                            "for information about why quantification failed. If this "
                            "problem persists, please report this issue on GitHub.");
            return 1;
        }

        bfs::path estFilePath = outputDirectory / "quant.sf";

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("geneMap")) {
            try {
                salmon::utils::generateGeneLevelEstimates(geneMapPath,
                                                            outputDirectory);
            } catch (std::exception& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
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
        std::cerr << "For usage information, " <<
            "try " << argv[0] << " quant --help-alignments\nExiting.\n";
        std::exit(1);
    }
    return 0;
}
