/**
>HEADER
    Copyright (c) 2013, 2014, 2015 Rob Patro rob.patro@cs.stonybrook.edu

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
#include <unordered_map>
#include <map>
#include <vector>
#include <unordered_set>
#include <iterator>
#include <mutex>
#include <thread>
#include <sstream>
#include <exception>
#include <random>
#include <queue>
#include <unordered_map>
#include "btree_map.h"
#include "btree_set.h"

// C++ string formatting library
#include "spdlog/details/format.h"

// C Includes for BWA
#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <cctype>

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "ksort.h"
#include "kvec.h"
#include "utils.h"
}

// Jellyfish 2 include
#include "jellyfish/mer_dna.hpp"
#include "jellyfish/stream_manager.hpp"
#include "jellyfish/whole_sequence_parser.hpp"

// Boost Includes
#include <boost/filesystem.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/dynamic_bitset/dynamic_bitset.hpp>
#include <boost/range/irange.hpp>
#include <boost/program_options.hpp>
#include <boost/lockfree/queue.hpp>
#include <boost/thread/thread.hpp>
#include <boost/range/iterator_range.hpp>

// TBB Includes
#include "tbb/concurrent_unordered_set.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_queue.h"
#include "tbb/parallel_for.h"
#include "tbb/parallel_for_each.h"
#include "tbb/parallel_reduce.h"
#include "tbb/blocked_range.h"
#include "tbb/task_scheduler_init.h"
#include "tbb/partitioner.h"

// logger includes
#include "spdlog/spdlog.h"

// Cereal includes
#include "cereal/types/vector.hpp"
#include "cereal/archives/binary.hpp"

#include "concurrentqueue.h"

// salmon / Salmon includes
#include "ClusterForest.hpp"
#include "SalmonMath.hpp"
#include "Transcript.hpp"
#include "LibraryFormat.hpp"
#include "SalmonUtils.hpp"
#include "ReadLibrary.hpp"
#include "SalmonConfig.hpp"
#include "IOUtils.hpp"
#include "SalmonIndex.hpp"

#include "BWAUtils.hpp"
#include "KmerIntervalMap.hpp"
#include "AlignmentGroup.hpp"
#include "PairSequenceParser.hpp"
#include "ForgettingMassCalculator.hpp"
#include "FragmentLengthDistribution.hpp"
#include "ReadExperiment.hpp"
#include "SalmonOpts.hpp"
#include "EquivalenceClassBuilder.hpp"
#include "CollapsedEMOptimizer.hpp"
#include "CollapsedGibbsSampler.hpp"
#include "RapMapUtils.hpp"
#include "HitManager.hpp"
#include "SASearcher.hpp"
#include "SACollector.hpp"
#include "TextBootstrapWriter.hpp"

/* This allows us to use CLASP for optimal MEM
 * chaining.  However, this seems to be neither
 * computationally efficient, nor significantly
 * better than the greedy chaining, so I'm temporarily
 * removing this un-necessary dependency.  If you
 * (other dev or future Rob) re-instate this in the future
 * remember to re-enable the CLASP fetch and build
 * steps in the CMakeLists.txt files
 *
 *#include "FragmentList.hpp"
 */

extern unsigned char nst_nt4_table[256];
char const* bwa_pg = "cha";

/****** QUASI MAPPING DECLARATIONS *********/
using MateStatus = rapmap::utils::MateStatus;
using QuasiAlignment = rapmap::utils::QuasiAlignment;
/****** QUASI MAPPING DECLARATIONS  *******/


/******* STUFF THAT IS STATIC IN BWAMEM THAT WE NEED HERE --- Just re-define it *************/
#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

typedef struct {
    bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

static smem_aux_t *smem_aux_init()
{
    smem_aux_t *a;
    a = static_cast<smem_aux_t*>(calloc(1, sizeof(smem_aux_t)));
    a->tmpv[0] = static_cast<bwtintv_v*>(calloc(1, sizeof(bwtintv_v)));
    a->tmpv[1] = static_cast<bwtintv_v*>(calloc(1, sizeof(bwtintv_v)));
    return a;
}

static void smem_aux_destroy(smem_aux_t *a)
{
    free(a->tmpv[0]->a); free(a->tmpv[0]);
    free(a->tmpv[1]->a); free(a->tmpv[1]);
    free(a->mem.a); free(a->mem1.a);
    free(a);
}

static void mem_collect_intv(const SalmonOpts& sopt, const mem_opt_t *opt, SalmonIndex* sidx, int len, const uint8_t *seq, smem_aux_t *a)
{
    const bwt_t* bwt = sidx->bwaIndex()->bwt;
    int i, k, x = 0, old_n;
    int start_width = (opt->flag & MEM_F_SELF_OVLP)? 2 : 1;
    int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
    a->mem.n = 0;

    // first pass: find all SMEMs
    if (sidx->hasAuxKmerIndex()) {
	std::cerr << "HERE1\n";
        KmerIntervalMap& auxIdx = sidx->auxIndex();
        uint32_t klen = auxIdx.k();
        while (x < len) {
            if (seq[x] < 4) {
                // Make sure there are at least k bases left
                if (len - x < klen) { x = len; continue; }
                // search for this key in the auxiliary index
                KmerKey kmer(const_cast<uint8_t*>(&(seq[x])), klen);
                auto it = auxIdx.find(kmer);
                // if we can't find it, move to the next key
                if (it == auxIdx.end()) { ++x; continue; }
                // otherwise, start the search using the initial interval @it->second from the hash
                int xb = x;
                x = bwautils::bwt_smem1_with_kmer(bwt, len, seq, x, start_width, it->second, &a->mem1, a->tmpv);
                for (i = 0; i < a->mem1.n; ++i) {
                    bwtintv_t *p = &a->mem1.a[i];
                    int slen = (uint32_t)p->info - (p->info>>32); // seed length
                    if (slen >= opt->min_seed_len)
                        kv_push(bwtintv_t, a->mem, *p);
                }
            } else ++x;
        }
    } else {
        while (x < len) {
            if (seq[x] < 4) {
                x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
                for (i = 0; i < a->mem1.n; ++i) {
                    bwtintv_t *p = &a->mem1.a[i];
                    int slen = (uint32_t)p->info - (p->info>>32); // seed length
                    if (slen >= opt->min_seed_len)
                        kv_push(bwtintv_t, a->mem, *p);
                }
            } else ++x;
        }
    }

    // For sensitive / extra-sensitive mode only
    if (sopt.sensitive or sopt.extraSeedPass) {
        // second pass: find MEMs inside a long SMEM
        old_n = a->mem.n;
        for (k = 0; k < old_n; ++k) {
            bwtintv_t *p = &a->mem.a[k];
            int start = p->info>>32, end = (int32_t)p->info;
            if (end - start < split_len || p->x[2] > opt->split_width) continue;

            //int idx = (start + end) >> 1;
            bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
            for (i = 0; i < a->mem1.n; ++i)
                if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
                    kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
        }
    }

    // For extra-sensitive mode only
    // third pass: LAST-like
    if (sopt.extraSeedPass and opt->max_mem_intv > 0) {
        x = 0;
        while (x < len) {
            if (seq[x] < 4) {
                if (1) {
                    bwtintv_t m;
                    x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
                    if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
                } else { // for now, we never come to this block which is slower
                    x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
                    for (i = 0; i < a->mem1.n; ++i)
                        kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
                }
            } else ++x;
        }
    }
    // sort
    // ks_introsort(mem_intv, a->mem.n, a->mem.a);
}


/******* END OF STUFF THAT IS STATIC IN BWAMEM THAT WE NEED HERE --- Just re-define it *************/

using paired_parser = pair_sequence_parser<char**>;
using stream_manager = jellyfish::stream_manager<std::vector<std::string>::const_iterator>;
using single_parser = jellyfish::whole_sequence_parser<stream_manager>;

using TranscriptID = uint32_t;
using TranscriptIDVector = std::vector<TranscriptID>;
using KmerIDMap = std::vector<TranscriptIDVector>;
using my_mer = jellyfish::mer_dna_ns::mer_base_static<uint64_t, 1>;

constexpr uint32_t miniBatchSize{5000};

class SMEMAlignment {
    public:
        SMEMAlignment() :
            transcriptID_(std::numeric_limits<TranscriptID>::max()),
            format_(LibraryFormat::formatFromID(0)),
            score_(0.0),
            hitPos_(0),
            fragLength_(0),
            logProb(salmon::math::LOG_0),
            logBias(salmon::math::LOG_0){}

        SMEMAlignment(TranscriptID transcriptIDIn, LibraryFormat format,
                  double scoreIn = 0.0,
                  int32_t hitPosIn = 0,
                  uint32_t fragLengthIn= 0,
                  double logProbIn = salmon::math::LOG_0) :
            transcriptID_(transcriptIDIn), format_(format), score_(scoreIn),
            hitPos_(hitPosIn), fragLength_(fragLengthIn), logProb(logProbIn) {}

        SMEMAlignment(const SMEMAlignment& o) = default;
        SMEMAlignment(SMEMAlignment&& o) = default;
        SMEMAlignment& operator=(SMEMAlignment& o) = default;
        SMEMAlignment& operator=(SMEMAlignment&& o) = default;


        inline TranscriptID transcriptID() const { return transcriptID_; }
        inline uint32_t fragLength() const { return fragLength_; }
        inline LibraryFormat libFormat() const { return format_; }
        inline double score() const { return score_; }
        inline int32_t hitPos() const { return hitPos_; }
        // inline double coverage() {  return static_cast<double>(kmerCount) / fragLength_; };
        uint32_t kmerCount;
        double logProb;
        double logBias;
        template <typename Archive>
        void save(Archive& archive) const {
            archive(transcriptID_, format_.formatID(), score_, hitPos_, fragLength_);
        }

        template <typename Archive>
        void load(Archive& archive) {
            uint8_t formatID;
            archive(transcriptID_, formatID, score_, hitPos_, fragLength_);
            format_ = LibraryFormat::formatFromID(formatID);
        }

    private:
        TranscriptID transcriptID_;
        LibraryFormat format_;
        double score_;
        int32_t hitPos_;
        uint32_t fragLength_;
};

template <typename AlnT>
using AlnGroupVec = std::vector<AlignmentGroup<AlnT>>;

template <typename AlnT>
using AlnGroupVecRange = boost::iterator_range<typename AlnGroupVec<AlnT>::iterator>;

#define __MOODYCAMEL__
#if defined(__MOODYCAMEL__)
 template <typename AlnT>
 using AlnGroupQueue = moodycamel::ConcurrentQueue<AlignmentGroup<AlnT>*>;
#else
 template <typename AlnT>
 using AlnGroupQueue = tbb::concurrent_queue<AlignmentGroup<AlnT>*>;
#endif

template <typename AlnT>
void processMiniBatch(
        ReadExperiment& readExp,
        ForgettingMassCalculator& fmCalc,
        uint64_t firstTimestepOfRound,
        ReadLibrary& readLib,
        const SalmonOpts& salmonOpts,
        //std::vector<AlignmentGroup<AlnT>*>& batchHits,
        AlnGroupVecRange<AlnT> batchHits,
        std::vector<Transcript>& transcripts,
        ClusterForest& clusterForest,
        FragmentLengthDistribution& fragLengthDist,
        std::atomic<uint64_t>& numAssignedFragments,
        std::default_random_engine& randEng,
        bool initialRound,
        bool& burnedIn
        ) {

    using salmon::math::LOG_0;
    using salmon::math::LOG_1;
    using salmon::math::LOG_ONEHALF;
    using salmon::math::logAdd;
    using salmon::math::logSub;

    const uint64_t numBurninFrags = salmonOpts.numBurninFrags;

    auto log = spdlog::get("jointLog");
    size_t numTranscripts{transcripts.size()};
    size_t localNumAssignedFragments{0};
    size_t priorNumAssignedFragments{numAssignedFragments};
    std::uniform_real_distribution<> uni(0.0, 1.0 + std::numeric_limits<double>::min());
    std::vector<uint64_t> libTypeCounts(LibraryFormat::maxLibTypeID() + 1);

    std::vector<FragmentStartPositionDistribution>& fragStartDists =
        readExp.fragmentStartPositionDistributions();
    auto& biasModel = readExp.sequenceBiasModel();

    bool updateCounts = initialRound;
    bool useReadCompat = salmonOpts.incompatPrior != salmon::math::LOG_1;
    bool useFSPD{!salmonOpts.noFragStartPosDist};
    bool useFragLengthDist{!salmonOpts.noFragLengthDist};
    bool useSequenceBiasModel{!salmonOpts.noSeqBiasModel};
    bool noFragLenFactor{salmonOpts.noFragLenFactor};

    const auto expectedLibraryFormat = readLib.format();
    uint64_t zeroProbFrags{0};

    //EQClass
    EquivalenceClassBuilder& eqBuilder = readExp.equivalenceClassBuilder();

    // Build reverse map from transcriptID => hit id
    using HitID = uint32_t;
    /* This isn't used anymore!!!
    btree::btree_map<TranscriptID, std::vector<SMEMAlignment*>> hitsForTranscript;
    size_t hitID{0};
    for (auto& hv : batchHits) {
        for (auto& tid : hv->alignments()) {
            hitsForTranscript[tid.transcriptID()].push_back(&tid);
        }
        ++hitID;
    }
    double clustTotal = std::log(batchHits.size()) + logForgettingMass;
    */

    double logForgettingMass{0.0};
    uint64_t currentMinibatchTimestep{0};
    fmCalc.getLogMassAndTimestep(logForgettingMass, currentMinibatchTimestep);

    double startingCumulativeMass = fmCalc.cumulativeLogMassAt(firstTimestepOfRound);
    // BEGIN: DOUBLY-COLLAPSED TESTING
    struct HitInfo {
        uint32_t numHits = 0;
        bool observed = false;
        double newUniqueMass = LOG_0;
    };

    int i{0};
    {
        // Iterate over each group of alignments (a group consists of all alignments reported
        // for a single read).  Distribute the read's mass to the transcripts
        // where it potentially aligns.
        for (auto& alnGroup : batchHits) {
            if (alnGroup.size() == 0) { continue; }

            // We start out with probability 0
            double sumOfAlignProbs{LOG_0};
            // Record whether or not this read is unique to a single transcript.
            bool transcriptUnique{true};

            auto firstTranscriptID = alnGroup.alignments().front().transcriptID();
            std::unordered_set<size_t> observedTranscripts;

            // EQCLASS
            std::vector<uint32_t> txpIDs;
            std::vector<double> auxProbs;
            double auxDenom = salmon::math::LOG_0;
	        size_t txpIDsHash{0};

            double avgLogBias = salmon::math::LOG_0;
            uint32_t numInGroup{0};
            // For each alignment of this read
            for (auto& aln : alnGroup.alignments()) {
                auto transcriptID = aln.transcriptID();
                auto& transcript = transcripts[transcriptID];
                transcriptUnique = transcriptUnique and (transcriptID == firstTranscriptID);

                double refLength = transcript.RefLength > 0 ? transcript.RefLength : 1.0;
                double coverage = aln.score();
                double logFragCov = (coverage > 0) ? std::log(coverage) : LOG_0;

                // The alignment probability is the product of a
                // transcript-level term (based on abundance and) an
                // alignment-level term.
                double logRefLength;
                if (salmonOpts.noEffectiveLengthCorrection or !burnedIn) {
                    logRefLength = std::log(transcript.RefLength);
                } else {
                    logRefLength = transcript.getCachedEffectiveLength();
                }

                double transcriptLogCount = transcript.mass(initialRound);

                if ( transcriptLogCount != LOG_0 ) {
                    double errLike = salmon::math::LOG_1;
                    if (burnedIn) {
                        // TODO: Make error model for smem-based quantification
                        //errLike = errMod.logLikelihood(aln, transcript);
                    }

                    double logFragProb = (useFragLengthDist) ?
                        ((aln.fragLength() > 0) ?
                         fragLengthDist.pmf(static_cast<size_t>(aln.fragLength())) :
                         LOG_1) :
                         LOG_1;

                    // TESTING
                    if (noFragLenFactor) { logFragProb = LOG_1; }


                    // TODO: Take the fragment length distribution into account
                    // for single-end fragments as in the alignment-based code below
                    /*
                    if (!salmonOpts.noFragLengthDist) {
                        if(aln->fragLen() == 0) {
                            if (aln->isLeft() and transcript.RefLength - aln->left() < fragLengthDist.maxVal()) {
                                logFragProb = fragLengthDist.cmf(transcript.RefLength - aln->left());
                            } else if (aln->isRight() and aln->right() < fragLengthDist.maxVal()) {
                                logFragProb = fragLengthDist.cmf(aln->right());
                            }
                        } else {
                            logFragProb = fragLengthDist.pmf(static_cast<size_t>(aln->fragLen()));
                        }
                    }
                    */

                    // The probability that the fragments align to the given strands in the
                    // given orientations.
                    double logAlignCompatProb = (useReadCompat) ?
                                                (salmon::utils::logAlignFormatProb(aln.libFormat(), expectedLibraryFormat, salmonOpts.incompatPrior)) :
                                                LOG_1;

                    // Allow for a non-uniform fragment start position distribution
                    double startPosProb = -logRefLength;
                    auto hitPos = aln.hitPos();
                    if (useFSPD and burnedIn and hitPos < refLength) {
                        auto& fragStartDist =
                            fragStartDists[transcript.lengthClassIndex()];
                        startPosProb = fragStartDist(hitPos, refLength, logRefLength);
                    }

                    double logBiasProb = salmon::math::LOG_1;
                    if (useSequenceBiasModel and burnedIn) {
                        double fragLength = aln.fragLength();
                        if (fragLength > 0) {
                            int32_t leftHitPos = hitPos;
                            int32_t rightHitPos = hitPos + fragLength;
                            logBiasProb = biasModel.biasFactor(transcript,
                                                               leftHitPos,
                                                               rightHitPos,
                                                               aln.libFormat());
                        } else {
                            logBiasProb = biasModel.biasFactor(transcript,
                                                               hitPos,
                                                               aln.libFormat());
                        }

                    }

                    // Increment the count of this type of read that we've seen
                    ++libTypeCounts[aln.libFormat().formatID()];

                    double auxProb = startPosProb + logFragProb + logFragCov +
                                     logAlignCompatProb + logBiasProb;

                    aln.logProb = transcriptLogCount + auxProb;

                    if (std::abs(aln.logProb) == LOG_0) { continue; }

                    if (useSequenceBiasModel and burnedIn) {
                        avgLogBias = salmon::math::logAdd(avgLogBias, logBiasProb);
                        numInGroup++;
                        aln.logBias = logBiasProb;
                    } else {
                        avgLogBias = salmon::math::logAdd(avgLogBias, logBiasProb);
                        numInGroup++;
                        aln.logBias = salmon::math::LOG_1;
                    }

                    sumOfAlignProbs = logAdd(sumOfAlignProbs, aln.logProb);

                    if (updateCounts and
                        observedTranscripts.find(transcriptID) == observedTranscripts.end()) {
                        transcripts[transcriptID].addTotalCount(1);
                        observedTranscripts.insert(transcriptID);
                    }
                    // EQCLASS
                    txpIDs.push_back(transcriptID);
                    auxProbs.push_back(auxProb);
                    auxDenom = salmon::math::logAdd(auxDenom, auxProb);
    	            boost::hash_combine(txpIDsHash, transcriptID);
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
            }

            if (numInGroup > 0){
                avgLogBias = avgLogBias - std::log(numInGroup);
            }

            // EQCLASS
            TranscriptGroup tg(txpIDs, txpIDsHash);
            double auxProbSum{0.0};
            for (auto& p : auxProbs) {
                p = std::exp(p - auxDenom);
                auxProbSum += p;
            }
            /*
            if (std::abs(auxProbSum - 1.0) > 0.01) {
                std::cerr << "weights had sum of " << auxProbSum
                          << " but it should be 1!!\n\n";
            }
            */
            eqBuilder.addGroup(std::move(tg), auxProbs);

            // normalize the hits
            for (auto& aln : alnGroup.alignments()) {
                if (std::abs(aln.logProb) == LOG_0) { continue; }
                // Normalize the log-probability of this alignment
                aln.logProb -= sumOfAlignProbs;
                // Get the transcript referenced in this alignment
                auto transcriptID = aln.transcriptID();
                auto& transcript = transcripts[transcriptID];

                // Add the new mass to this transcript
                double newMass = logForgettingMass + aln.logProb;
                transcript.addMass( newMass );
                //transcript.addBias( aln.logBias );

                double r = uni(randEng);
                if (!burnedIn and r < std::exp(aln.logProb)) {
                    //errMod.update(aln, transcript, aln.logProb, logForgettingMass);
                    double fragLength = aln.fragLength();
                    if (useFragLengthDist and fragLength > 0.0) {
                        //if (aln.fragType() == ReadType::PAIRED_END) {
                        fragLengthDist.addVal(fragLength, logForgettingMass);
                    }
                    if (useFSPD) {
                        auto hitPos = aln.hitPos();
                        auto& fragStartDist =
                            fragStartDists[transcript.lengthClassIndex()];
                        fragStartDist.addVal(hitPos,
                                             transcript.RefLength,
                                             logForgettingMass);
                    }
                    if (useSequenceBiasModel) {
                        if (fragLength > 0.0) {
                            int32_t leftPos = aln.hitPos();
                            int32_t rightPos = leftPos + fragLength;
                            biasModel.update(transcript, leftPos, rightPos,
                                             aln.libFormat(), logForgettingMass, LOG_1);
                        } else {
                            int32_t hitPos = aln.hitPos();
                            biasModel.update(transcript, hitPos,
                                             aln.libFormat(),
                                             logForgettingMass, LOG_1);
                        }
                    }
                }
            } // end normalize

            //double avgBias = std::exp(avgLogBias);
            // update the single target transcript
            if (transcriptUnique) {
                if (updateCounts) {
                    transcripts[firstTranscriptID].addUniqueCount(1);
                }
                clusterForest.updateCluster(
                        firstTranscriptID,
                        1.0,
                        logForgettingMass, updateCounts);
            } else { // or the appropriate clusters
                clusterForest.mergeClusters<AlnT>(alnGroup.alignments().begin(), alnGroup.alignments().end());
                clusterForest.updateCluster(
                        alnGroup.alignments().front().transcriptID(),
                        1.0,
                        logForgettingMass, updateCounts);
            }

            } // end read group
        }// end timer

        double individualTotal = LOG_0;
        {
            /*
            // M-step
            double totalMass{0.0};
            for (auto kv = hitsForTranscript.begin(); kv != hitsForTranscript.end(); ++kv) {
                auto transcriptID = kv->first;
                // The target must be a valid transcript
                if (transcriptID >= numTranscripts or transcriptID < 0) {std::cerr << "index " << transcriptID << " out of bounds\n"; }

                auto& transcript = transcripts[transcriptID];

                // The prior probability
                double hitMass{LOG_0};

                // The set of alignments that match transcriptID
                auto& hits = kv->second;
                std::for_each(hits.begin(), hits.end(), [&](SMEMAlignment* aln) -> void {
                        if (!std::isfinite(aln->logProb)) { std::cerr << "hitMass = " << aln->logProb << "\n"; }
                        hitMass = logAdd(hitMass, aln->logProb);
                        });

                double updateMass = logForgettingMass + hitMass;
                individualTotal = logAdd(individualTotal, updateMass);
                totalMass = logAdd(totalMass, updateMass);
                transcript.addMass(updateMass);
            } // end for
            */
        } // end timer

        if (zeroProbFrags > 0) {
            log->warn("Minibatch contained {} "
                      "0 probability fragments", zeroProbFrags);
        }

        numAssignedFragments += localNumAssignedFragments;
        if (numAssignedFragments >= numBurninFrags and !burnedIn) {
            burnedIn = true;
            for (auto& t : transcripts) {  t.updateEffectiveLength(fragLengthDist); }
            if (useFSPD) {
                // update all of the fragment start position
                // distributions
                for (auto& fspd : fragStartDists) {
                    fspd.update();
                }
            }
        }
        if (initialRound) {
            readLib.updateLibTypeCounts(libTypeCounts);
        }
}

uint32_t basesCovered(std::vector<uint32_t>& kmerHits) {
    std::sort(kmerHits.begin(), kmerHits.end());
    uint32_t covered{0};
    uint32_t lastHit{0};
    uint32_t kl{20};
    for (auto h : kmerHits) {
        covered += std::min(h - lastHit, kl);
        lastHit = h;
    }
    return covered;
}

uint32_t basesCovered(std::vector<uint32_t>& posLeft, std::vector<uint32_t>& posRight) {
    return basesCovered(posLeft) + basesCovered(posRight);
}

class KmerVote {
    public:
        KmerVote(int32_t vp, uint32_t rp, uint32_t vl) : votePos(vp), readPos(rp), voteLen(vl) {}
        int32_t votePos{0};
        uint32_t readPos{0};
        uint32_t voteLen{0};
        /*
        std::string str(){
            return "<" + votePos  + ", "  + readPos  + ", "  + voteLen + ">";
        }
        */
};
class MatchFragment {
    public:
        MatchFragment(uint32_t refStart_, uint32_t queryStart_, uint32_t length_) :
            refStart(refStart_), queryStart(queryStart_), length(length_) {}

        uint32_t refStart, queryStart, length;
        uint32_t weight;
        double score;
};

bool precedes(const MatchFragment& a, const MatchFragment& b) {
    return (a.refStart + a.length) < b.refStart and
           (a.queryStart + a.length) < b.queryStart;
}


class TranscriptHitList {
    public:
        int32_t bestHitPos{0};
        uint32_t bestHitCount{0};
        double bestHitScore{0.0};

        std::vector<KmerVote> votes;
        std::vector<KmerVote> rcVotes;

        uint32_t targetID;
        uint32_t fwdCov{0};
        uint32_t revCov{0};

        bool isForward_{true};

        void addFragMatch(uint32_t tpos, uint32_t readPos, uint32_t voteLen) {
            int32_t votePos = static_cast<int32_t>(tpos) - static_cast<int32_t>(readPos);
            votes.emplace_back(votePos, readPos, voteLen);
            fwdCov += voteLen;
        }

        void addFragMatchRC(uint32_t tpos, uint32_t readPos, uint32_t voteLen, uint32_t readLen) {
            //int32_t votePos = static_cast<int32_t>(tpos) - (readPos) + voteLen;
            int32_t votePos = static_cast<int32_t>(tpos) - (readLen - readPos);
            rcVotes.emplace_back(votePos, readPos, voteLen);
            revCov += voteLen;
        }

        uint32_t totalNumHits() { return std::max(votes.size(), rcVotes.size()); }

        bool computeBestLocFast_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                                 std::string& read, bool isRC,
                                 int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {
            bool updatedMaxScore{true};
            if (sVotes.size() == 0) { return updatedMaxScore; }
            uint32_t readLen = read.length();
            uint32_t votePos = sVotes.front().votePos;

            uint32_t cov = isRC ? revCov : fwdCov;
            if (cov > maxClusterCount) {
                maxClusterCount = cov;
                maxClusterPos = votePos;
                maxClusterScore = maxClusterCount / static_cast<double>(readLen);
                updatedMaxScore = true;
            }
            return updatedMaxScore;

        }

        bool computeBestLoc_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                             std::string& read, bool isRC,
                             int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {
            // Did we update the highest-scoring cluster? This will be set to
            // true iff we have a cluster of a higher score than the score
            // currently given in maxClusterCount.
            bool updatedMaxScore{false};

            if (sVotes.size() == 0) { return updatedMaxScore; }

            struct VoteInfo {
                uint32_t coverage = 0;
                int32_t rightmostBase = 0;
            };

            uint32_t readLen = read.length();

            boost::container::flat_map<uint32_t, VoteInfo> hitMap;
            int32_t currClust{static_cast<int32_t>(sVotes.front().votePos)};

            for (size_t j = 0; j < sVotes.size(); ++j) {

                int32_t votePos = sVotes[j].votePos;
                uint32_t readPos = sVotes[j].readPos;
                uint32_t voteLen = sVotes[j].voteLen;

                if (votePos >= currClust) {
                    if (votePos - currClust > 10) {
                        currClust = votePos;
                    }
                    auto& hmEntry = hitMap[currClust];

                    hmEntry.coverage += std::min(voteLen, (votePos + readPos + voteLen) - hmEntry.rightmostBase);
                    hmEntry.rightmostBase = votePos + readPos + voteLen;
                } else if (votePos < currClust) {
                    std::cerr << "Should not have votePos = " << votePos << " <  currClust = " << currClust << "\n";
                    std::exit(1);
                }

                if (hitMap[currClust].coverage > maxClusterCount) {
                    maxClusterCount = hitMap[currClust].coverage;
                    maxClusterPos = currClust;
                    maxClusterScore = maxClusterCount / static_cast<double>(readLen);
                    updatedMaxScore = true;
                }

            }
            return updatedMaxScore;
        }

        bool computeBestLoc2_(std::vector<KmerVote>& sVotes, uint32_t tlen,
                              int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {

            bool updatedMaxScore{false};

            if (sVotes.size() == 0) { return updatedMaxScore; }

            double weights[] = { 1.0, 0.983471453822, 0.935506985032,
                0.860707976425, 0.765928338365, 0.6592406302, 0.548811636094,
                0.441902209585, 0.344153786865, 0.259240260646,
                0.188875602838};

            uint32_t maxGap = 4;
            uint32_t leftmost = (sVotes.front().votePos > maxGap) ? (sVotes.front().votePos - maxGap) : 0;
            uint32_t rightmost = std::min(sVotes.back().votePos + maxGap, tlen);

            uint32_t span = (rightmost - leftmost);
            std::vector<double> probAln(span, 0.0);
            double kwidth = 1.0 / (2.0 * maxGap);

            size_t nvotes = sVotes.size();
            for (size_t j = 0; j < nvotes; ++j) {
                uint32_t votePos = sVotes[j].votePos;
                uint32_t voteLen = sVotes[j].voteLen;

                auto x = j + 1;
                while (x < nvotes and sVotes[x].votePos == votePos) {
                    voteLen += sVotes[x].voteLen;
                    j += 1;
                    x += 1;
                }


                uint32_t dist{0};
                size_t start = (votePos >= maxGap) ? (votePos - maxGap - leftmost) : (votePos - leftmost);
                size_t mid = votePos - leftmost;
                size_t end = std::min(votePos + maxGap - leftmost, rightmost - leftmost);
                for (size_t k = start; k < end; k += 1) {
                    dist = (mid > k) ? mid - k : k - mid;
                    probAln[k] += weights[dist] * voteLen;
                    if (probAln[k] > maxClusterScore) {
                        maxClusterScore = probAln[k];
                        maxClusterPos = k + leftmost;
                        updatedMaxScore = true;
                    }
                }
            }

            return updatedMaxScore;
        }


        inline uint32_t numSampledHits_(Transcript& transcript, std::string& readIn,
                                        int32_t votePos, int32_t posInRead, int32_t voteLen, bool isRC, uint32_t numTries) {


            // The read starts at this position in the transcript (may be negative!)
            int32_t readStart = votePos;
            // The (uncorrected) length of the read
            int32_t readLen = readIn.length();
            // Pointer to the sequence of the read
            const char* read = readIn.c_str();
            // Don't mess around with unsigned arithmetic here
            int32_t tlen = transcript.RefLength;

            // If the read starts before the first base of the transcript,
            // trim off the initial overhang  and correct the other variables
            if (readStart < 0) {
                if (isRC) {
                    uint32_t correction = -readStart;
                    //std::cerr << "readLen = " << readLen << ", posInRead = " << posInRead << ", voteLen = " << voteLen << ", correction = " << correction << "\n";
                    //std::cerr << "tlen = " << tlen << ", votePos = " << votePos << "\n";
                    read += correction;
                    readLen -= correction;
                    posInRead -= correction;
                    readStart = 0;
                } else {
                    uint32_t correction = -readStart;
                    read += correction;
                    readLen -= correction;
                    posInRead -= correction;
                    readStart = 0;
                }
            }
            // If the read hangs off the end of the transcript,
            // shorten its effective length.
            if (readStart + readLen >= tlen) {
                if (isRC) {
                    uint32_t correction = (readStart + readLen) - transcript.RefLength + 1;
                    //std::cerr << "Trimming RC hit: correction = " << correction << "\n";
                    //std::cerr << "untrimmed read : "  << read << "\n";
                    read += correction;
                    readLen -= correction;
                    if (voteLen > readLen) { voteLen = readLen; }
                    posInRead = 0;
                } else {
                    readLen = tlen - (readStart + 1);
                    voteLen = std::max(voteLen, readLen - (posInRead + voteLen));
                }
            }
            // Finally, clip any reverse complement reads starting at 0
            if (isRC) {

                if (voteLen > readStart) {
                    readLen -= (readLen - (posInRead + voteLen));
                }

            }

            // If the read is too short, it's not useful
            if (readLen <= 15) { return 0; }
            // The step between sample centers (given the number of samples we're going to take)
            double step = (readLen - 1) / static_cast<double>(numTries-1);
            // The strand of the transcript from which we'll extract sequence
            auto dir = (isRC) ? salmon::stringtools::strand::reverse :
                                salmon::stringtools::strand::forward;

            bool superVerbose{false};

            if (superVerbose) {
                std::stringstream ss;
                ss << "Supposed hit " << (isRC ? "RC" : "") << "\n";
                ss << "info: votePos = " << votePos << ", posInRead = " << posInRead
                    << ", voteLen = " << voteLen << ", readLen = " << readLen
                    << ", tran len = " << tlen << ", step = " << step << "\n";
                if (readStart + readLen > tlen ) {
                    ss << "ERROR!!!\n";
                    std::cerr << "[[" << ss.str() << "]]";
                    std::exit(1);
                }
                ss << "Transcript name = " << transcript.RefName << "\n";
                ss << "T : ";
                try {
                    for ( size_t j = 0; j < readLen; ++j) {
                        if (isRC) {
                            if (j == posInRead) {
                                char red[] = "\x1b[30m";
                                red[3] = '0' + static_cast<char>(fmt::RED);
                                ss << red;
                            }

                            if (j == posInRead + voteLen) {
                                const char RESET_COLOR[] = "\x1b[0m";
                                ss << RESET_COLOR;
                            }
                            ss << transcript.charBaseAt(readStart+readLen-j,dir);
                        } else {
                            if (j == posInRead ) {
                                char red[] = "\x1b[30m";
                                red[3] = '0' + static_cast<char>(fmt::RED);
                                ss << red;
                            }

                            if (j == posInRead + voteLen) {
                                const char RESET_COLOR[] = "\x1b[0m";
                                ss << RESET_COLOR;
                            }

                            ss << transcript.charBaseAt(readStart+j);
                        }
                    }
                    ss << "\n";
                    char red[] = "\x1b[30m";
                    red[3] = '0' + static_cast<char>(fmt::RED);
                    const char RESET_COLOR[] = "\x1b[0m";

                    ss << "R : " << std::string(read, posInRead) << red << std::string(read + posInRead, voteLen) << RESET_COLOR;
                    if (readLen > posInRead + voteLen) { ss << std::string(read + posInRead + voteLen); }
                    ss << "\n\n";
                } catch (std::exception& e) {
                    std::cerr << "EXCEPTION !!!!!! " << e.what() << "\n";
                }
                std::cerr << ss.str() << "\n";
                ss.clear();
            }

            // The index of the current sample within the read
            int32_t readIndex = 0;

            // The number of loci in the subvotes and their
            // offset patternns
            size_t lpos = 3;
            int leftPattern[] = {-4, -2, 0};
            int rightPattern[] = {0, 2, 4};
            int centerPattern[] = {-4, 0, 4};

            // The number of subvote hits we've had
            uint32_t numHits = 0;
            // Take the samples
            for (size_t i  = 0; i < numTries; ++i) {
                // The sample will be centered around this point
                readIndex = static_cast<uint32_t>(std::round(readStart + i * step)) - readStart;

                // The number of successful sub-ovtes we have
                uint32_t subHit = 0;
                // Select the center sub-vote pattern, unless we're near the end of a read
                int* pattern = &centerPattern[0];
                if (readIndex + pattern[0] < 0) {
                    pattern = &rightPattern[0];
                } else if (readIndex + pattern[lpos-1] >= readLen) {
                    pattern = &leftPattern[0];
                }

                // collect the subvotes
                for (size_t j = 0; j < lpos; ++j) {
                    // the pattern offset
                    int offset = pattern[j];
                    // and sample position it implies within the read
                    int readPos = readIndex + offset;

                    if (readStart + readPos >= tlen) {
                        std::cerr  << "offset = " << offset << ", readPos = " << readPos << ", readStart = " << readStart << ", readStart + readPos = " << readStart + readPos << ", tlen = " << transcript.RefLength << "\n";
                    }

                    subHit += (isRC) ?
                        (transcript.charBaseAt(readStart + readLen - readPos, dir) == salmon::stringtools::charCanon[read[readPos]]) :
                        (transcript.charBaseAt(readStart + readPos               ) == salmon::stringtools::charCanon[read[readPos]]);
                }
                // if the entire subvote was successful, this is a hit
                numHits += (subHit == lpos);
            }
            // return the number of hits we had
            return numHits;
        }



        bool computeBestLoc3_(std::vector<KmerVote>& sVotes, Transcript& transcript,
                              std::string& read, bool isRC,
                              int32_t& maxClusterPos, uint32_t& maxClusterCount, double& maxClusterScore) {

            bool updatedMaxScore{false};

            if (sVotes.size() == 0) { return updatedMaxScore; }

            struct LocHitCount {
                int32_t loc;
                uint32_t nhits;
            };

            uint32_t numSamp = 15;
            std::vector<LocHitCount> hitCounts;
            size_t nvotes = sVotes.size();
            int32_t prevPos = -std::numeric_limits<int32_t>::max();
            for (size_t j = 0; j < nvotes; ++j) {
                int32_t votePos = sVotes[j].votePos;
                int32_t posInRead = sVotes[j].readPos;
                int32_t voteLen = sVotes[j].voteLen;
                if (prevPos == votePos) { continue; }
                auto numHits = numSampledHits_(transcript, read, votePos, posInRead, voteLen, isRC, numSamp);
                hitCounts.push_back({votePos, numHits});
                prevPos = votePos;
            }

            uint32_t maxGap = 8;
            uint32_t hitIdx = 0;
            uint32_t accumHits = 0;
            int32_t hitLoc = hitCounts[hitIdx].loc;
            while (hitIdx < hitCounts.size()) {
                uint32_t idx2 = hitIdx;
                while (idx2 < hitCounts.size() and std::abs(hitCounts[idx2].loc - hitLoc) <= maxGap) {
                    accumHits += hitCounts[idx2].nhits;
                    ++idx2;
                }

                double score = static_cast<double>(accumHits) / numSamp;
                if (score > maxClusterScore) {
                    maxClusterCount = accumHits;
                    maxClusterScore = score;
                    maxClusterPos = hitCounts[hitIdx].loc;
                    updatedMaxScore = true;
                }
                accumHits = 0;
                ++hitIdx;
                hitLoc = hitCounts[hitIdx].loc;
            }

            return updatedMaxScore;
        }


        bool computeBestChain(Transcript& transcript, std::string& read) {
            std::sort(votes.begin(), votes.end(),
                    [](const KmerVote& v1, const KmerVote& v2) -> bool {
                        if (v1.votePos == v2.votePos) {
                            return v1.readPos < v2.readPos;
                        }
                        return v1.votePos < v2.votePos;
                    });

            std::sort(rcVotes.begin(), rcVotes.end(),
                    [](const KmerVote& v1, const KmerVote& v2) -> bool {
                        if (v1.votePos == v2.votePos) {
                            return v1.readPos < v2.readPos;
                        }
                        return v1.votePos < v2.votePos;
                    });

            int32_t maxClusterPos{0};
            uint32_t maxClusterCount{0};
            double maxClusterScore{0.0};

            // we don't need the return value from the first call
            static_cast<void>(computeBestLoc_(votes, transcript, read, false, maxClusterPos, maxClusterCount, maxClusterScore));
            bool revIsBest = computeBestLoc_(rcVotes, transcript, read, true, maxClusterPos, maxClusterCount, maxClusterScore);
            isForward_ = not revIsBest;

            bestHitPos = maxClusterPos;
            bestHitCount = maxClusterCount;
            bestHitScore = maxClusterScore;
            return true;
        }

        bool isForward() { return isForward_; }

};

template <typename CoverageCalculator>
inline void collectHitsForRead(SalmonIndex* sidx, const bwtintv_v* a, smem_aux_t* auxHits,
                        mem_opt_t* memOptions, const SalmonOpts& salmonOpts, const uint8_t* read, uint32_t readLen,
                        std::vector<CoverageCalculator>& hits) {
                        //std::unordered_map<uint64_t, CoverageCalculator>& hits) {

    bwaidx_t* idx = sidx->bwaIndex();
    mem_collect_intv(salmonOpts, memOptions, sidx, readLen, read, auxHits);

    // For each MEM
    int firstSeedLen{-1};
    for (int i = 0; i < auxHits->mem.n; ++i ) {
        // A pointer to the interval of the MEMs occurences
        bwtintv_t* p = &auxHits->mem.a[i];
        // The start and end positions in the query string (i.e. read) of the MEM
        int qstart = p->info>>32;
        uint32_t qend = static_cast<uint32_t>(p->info);
        int step, count, slen = (qend - qstart); // seed length

        /*
        if (firstSeedLen > -1) {
            if (slen < firstSeedLen) { return; }
        } else {
            firstSeedLen = slen;
        }
        */

        int64_t k;
        step = p->x[2] > memOptions->max_occ? p->x[2] / memOptions->max_occ : 1;
        // For every occurrence of the MEM
        for (k = count = 0; k < p->x[2] && count < memOptions->max_occ; k += step, ++count) {
            bwtint_t pos;
            bwtint_t startPos, endPos;
            int len, isRev, isRevStart, isRevEnd, refID, refIDStart, refIDEnd;
            int queryStart = qstart;
            len = slen;
            uint32_t rlen = readLen;

            // Get the position in the reference index of this MEM occurrence
            int64_t refStart = bwt_sa(idx->bwt, p->x[0] + k);

            pos = startPos = bns_depos(idx->bns, refStart, &isRevStart);
            endPos = bns_depos(idx->bns, refStart + slen - 1, &isRevEnd);
            // If we span the forward/reverse boundary, discard the hit
            if (isRevStart != isRevEnd) {
                continue;
            }
            // Otherwise, isRevStart = isRevEnd so just assign isRev = isRevStart
            isRev = isRevStart;

            // If the hit is reversed --- swap the start and end
            if (isRev) {
                if (endPos > startPos) {
                    salmonOpts.jointLog->warn("Hit is supposedly reversed, "
                                              "but startPos = {} < endPos = {}",
                                              startPos, endPos);
                }
                auto temp = startPos;
                startPos = endPos;
                endPos = temp;
            }
            // Get the ID of the reference sequence in which it occurs
            refID = refIDStart = bns_pos2rid(idx->bns, startPos);
            refIDEnd = bns_pos2rid(idx->bns, endPos);

            if (refID < 0) { continue; } // bridging multiple reference sequences or the forward-reverse boundary;

            auto tlen = idx->bns->anns[refID].len;

            // The refence sequence-relative (e.g. transcript-relative) position of the MEM
            long hitLoc = static_cast<long>(isRev ? endPos : startPos) - idx->bns->anns[refID].offset;

            if ((refIDStart != refIDEnd)) {
                // If a seed spans two transcripts

                // If we're not considering splitting such seeds, then
                // just discard this seed and continue.
                if (not salmonOpts.splitSpanningSeeds) { continue; }

                //std::cerr << "Seed spans two transcripts! --- attempting to split: \n";
                if (!isRev) {
                    // If it's going forward, we have a situation like this
                    // packed transcripts: t1 ===========|t2|==========>
                    // hit:                          |==========>

                    // length of hit in t1
                    auto len1 = tlen - hitLoc;
                    // length of hit in t2
                    auto len2 = slen - len1;
                    if (std::max(len1, len2) < memOptions->min_seed_len) { continue; }

                    /** Keeping this here for now in case I need to debug splitting seeds again
                    std::cerr << "\t hit is in the forward direction: ";
                    std::cerr << "t1 part has length " << len1 << ", t2 part has length " << len2 << "\n";
                    */

                    // If the part in t1 is larger then just cut off the rest
                    if (len1 >= len2) {
                        slen = len1;
                        int32_t votePos = static_cast<int32_t>(hitLoc) - queryStart;
                        //std::cerr << "\t\t t1 (of length " << tlen << ") has larger hit --- new hit length = " << len1 << "; starts at pos " << queryStart << " in the read (votePos will be " << votePos << ")\n";
                    } else {
                        // Otherwise, make the hit be in t2.
                        // Because the hit spans the boundary where t2 begins,
                        // the new seed begins matching at position 0 of
                        // transcript t2
                        hitLoc = 0;
                        slen = len2;
                        // The seed originally started at position q, now it starts  len1 characters to the  right of that
                        queryStart += len1;
                        refID = refIDEnd;
                        int32_t votePos = static_cast<int32_t>(hitLoc) - queryStart;
                        tlen = idx->bns->anns[refID].len;
                        //std::cerr << "\t\t t2 (of length " << tlen << ") has larger hit --- new hit length = " << len2 << "; starts at pos " << queryStart << " in the read (votePos will be " << votePos << ")\n";
                    }
                } else {

                    // If it's going in the reverse direction, we have a situation like this
                    // packed transcripts: t1 <===========|t2|<==========
                    // hit:                          X======Y>======Z>
                    // Which means we have
                    // packed transcripts: t1 <===========|t2|<==========
                    // hit:                          <Z=====Y<======X
                    // length of hit in t1

                    auto len2 = endPos - idx->bns->anns[refIDEnd].offset;
                    auto len1 = slen - len2;
                    if (std::max(len1, len2) < memOptions->min_seed_len) { continue; }

                    /** Keeping this here for now in case I need to debug splitting seeds again
                    std::cerr << "\t hit is in the reverse direction: ";
                    std::cerr << "\n\n";
                    std::cerr << "startPos = " << startPos << ", endPos = " << endPos << ", offset[refIDStart] = "
                              <<  idx->bns->anns[refIDStart].offset << ", offset[refIDEnd] = " << idx->bns->anns[refIDEnd].offset << "\n";
                    std::cerr << "\n\n";
                    std::cerr << "t1 part has length " << len1 << ", t2 part has length " << len2 << "\n\n";
                    */

                    if (len1 >= len2) {
                        slen = len1;
                        hitLoc = tlen - len2;
                        queryStart += len2;
                        rlen -= len2;
                        int32_t votePos = static_cast<int32_t>(hitLoc) - (rlen - queryStart);
                        //std::cerr << "\t\t t1 (hitLoc: " << hitLoc << ") (of length " << tlen << ") has larger hit --- new hit length = " << len1 << "; starts at pos " << queryStart << " in the read (votePos will be " << votePos << ")\n";
                    } else {
                        slen = len2;
                        refID = bns_pos2rid(idx->bns, endPos);
                        tlen = idx->bns->anns[refID].len;
                        hitLoc = len2;
                        rlen = hitLoc + queryStart;
                        int32_t votePos = static_cast<int32_t>(hitLoc) - (rlen - queryStart);
                        //std::cerr << "\t\t t2 (of length " << tlen << ") (hitLoc: " << hitLoc << ") has larger hit --- new hit length = " << len2 << "; starts at pos " << queryStart << " in the read (votePos will be " << votePos << ")\n";
                    }
                }

            }

            auto hitIt = std::find_if(hits.begin(), hits.end(), [refID](CoverageCalculator& c) -> bool { return c.targetID == refID; });
            if (isRev) {
                if (hitIt == hits.end()) {
                    CoverageCalculator hit;
                    hit.targetID = refID;
                    hit.addFragMatchRC(hitLoc, queryStart, slen, rlen);
                    hits.emplace_back(hit);
                } else {
                    hitIt->addFragMatchRC(hitLoc, queryStart , slen, rlen);
                    //hits[refID].addFragMatchRC(hitLoc, queryStart , slen, rlen);
                }
            } else {
                if (hitIt == hits.end()) {
                    CoverageCalculator hit;
                    hit.targetID = refID;
                    hit.addFragMatch(hitLoc, queryStart, slen);
                    hits.emplace_back(hit);
                } else {
                    hitIt->addFragMatch(hitLoc, queryStart , slen);
                    //hits[refID].addFragMatch(hitLoc, queryStart, slen);
                }
            }
        } // for k
    }
}

inline bool consistentNames(header_sequence_qual& r) {
    return true;
}

bool consistentNames(std::pair<header_sequence_qual, header_sequence_qual>& rp) {
        auto l1 = rp.first.header.length();
        auto l2 = rp.second.header.length();
        char* sptr = static_cast<char*>(memchr(&rp.first.header[0], ' ', l1));

        bool compat = false;
        // If we didn't find a space in the name of read1
        if (sptr == NULL) {
            if (l1 > 1) {
                compat = (l1 == l2);
                compat = compat and (memcmp(&rp.first.header[0], &rp.second.header[0], l1-1) == 0);
                compat = compat and ((rp.first.header[l1-1] == '1' and rp.second.header[l2-1] == '2')
                                or   (rp.first.header[l1-1] == rp.second.header[l2-1]));
            } else {
                compat = (l1 == l2);
                compat = compat and (rp.first.header[0] == rp.second.header[0]);
            }
        } else {
            size_t offset = sptr - (&rp.first.header[0]);

            // If read2 matches read1 up to and including the space
            if (offset + 1 < l2) {
                compat = memcmp(&rp.first.header[0], &rp.second.header[0], offset) == 0;
                // and after the space, read1 and read2 have an identical character or
                // read1 has a '1' and read2 has a '2', then this is a consistent pair.
                compat = compat and ((rp.first.header[offset+1] == rp.second.header[offset+1])
                                or   (rp.first.header[offset+1] == '1' and rp.second.header[offset+1] == '2'));
            } else {
                compat = false;
            }
        }
        return compat;
}

/**
 *  Returns true if the @hit is within @cutoff bases of the end of
 *  transcript @txp and false otherwise.
 */
template <typename CoverageCalculator>
inline bool nearEndOfTranscript(
            CoverageCalculator& hit,
            Transcript& txp,
            int32_t cutoff=std::numeric_limits<int32_t>::max()) {
	// check if hit appears close to the end of the given transcript
    bool isForward = hit.isForward();
	int32_t hitPos = static_cast<int32_t>(hit.bestHitPos);
    return (hitPos <= cutoff or std::abs(static_cast<int32_t>(txp.RefLength) - hitPos) <= cutoff);
}

template <typename CoverageCalculator>
inline void getHitsForFragment(std::pair<header_sequence_qual, header_sequence_qual>& frag,
                        SalmonIndex* sidx,
                        smem_i *itr,
                        const bwtintv_v *a,
                        smem_aux_t* auxHits,
                        mem_opt_t* memOptions,
                        ReadExperiment& readExp,
                        const SalmonOpts& salmonOpts,
                        double coverageThresh,
                        uint64_t& upperBoundHits,
                        AlignmentGroup<SMEMAlignment>& hitList,
                        uint64_t& hitListCount,
                        std::vector<Transcript>& transcripts) {

    //std::unordered_map<uint64_t, CoverageCalculator> leftHits;
    //std::unordered_map<uint64_t, CoverageCalculator> rightHits;

    std::vector<CoverageCalculator> leftHits;
    std::vector<CoverageCalculator> rightHits;


    uint32_t leftReadLength{0};
    uint32_t rightReadLength{0};

    auto& eqBuilder = readExp.equivalenceClassBuilder();
    bool allowOrphans{salmonOpts.allowOrphans};

    /**
    * As soon as we can decide on an acceptable way to validate read names,
    * we'll inform the user and quit if we see something inconsistent.  However,
    * we first need a reasonable way to verify potential naming formats from
    * many different sources.
    */
    /*
    if (!consistentNames(frag)) {
        fmt::MemoryWriter errstream;

        errstream << "Inconsistent paired-end reads!\n";
        errstream << "mate1 : " << frag.first.header << "\n";
        errstream << "mate2 : " << frag.second.header << "\n";
        errstream << "Paired-end reads should appear consistently in their respective files.\n";
        errstream << "Please fix the paire-end input before quantifying with salmon; exiting.\n";

        std::cerr << errstream.str();
        std::exit(-1);
    }
    */

    //---------- End 1 ----------------------//
    {
        std::string readStr   = frag.first.seq;
        uint32_t readLen      = readStr.size();

        leftReadLength = readLen;

        for (int p = 0; p < readLen; ++p) {
            readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
        }

        collectHitsForRead(sidx, a, auxHits,
                            memOptions,
                            salmonOpts,
                            reinterpret_cast<const uint8_t*>(readStr.c_str()),
                            readLen,
                            leftHits);
    }

    //---------- End 2 ----------------------//
    {
        std::string readStr   = frag.second.seq;
        uint32_t readLen      = readStr.size();

        rightReadLength = readLen;

        for (int p = 0; p < readLen; ++p) {
            readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
        }

        collectHitsForRead(sidx, a, auxHits,
                            memOptions,
                            salmonOpts,
                            reinterpret_cast<const uint8_t*>(readStr.c_str()),
                            readLen,
                            rightHits);
     } // end right

    size_t numTrivialHits = (leftHits.size() + rightHits.size() > 0) ? 1 : 0;
    upperBoundHits += (leftHits.size() + rightHits.size() > 0) ? 1 : 0;
    size_t readHits{0};
    auto& alnList = hitList.alignments();
    hitList.isUniquelyMapped() = true;
    alnList.clear();
    // nothing more to do
    if (numTrivialHits == 0) { return; }


    double cutoffLeft{ coverageThresh };//* leftReadLength};
    double cutoffRight{ coverageThresh };//* rightReadLength};

    uint64_t leftHitCount{0};

    // Fraction of the optimal coverage that a lightweight alignment
    // must obtain in order to be retained.
    float fOpt{0.95};

    // First, see if there are transcripts where both ends of the
    // fragments map
    auto& minHitList = (leftHits.size() < rightHits.size()) ? leftHits : rightHits;
    auto& maxHitList = (leftHits.size() < rightHits.size()) ? rightHits : leftHits;

    struct JointHitPtr {
        uint32_t transcriptID;
        size_t leftIndex;
        size_t rightIndex;
    };

    std::vector<JointHitPtr> jointHits; // haha (variable name)!
    jointHits.reserve(minHitList.size());

    // vector-based code
    // Sort the left and right hits
    std::sort(leftHits.begin(), leftHits.end(),
              [](CoverageCalculator& c1, CoverageCalculator& c2) -> bool {
                return c1.targetID < c2.targetID;
               });
    std::sort(rightHits.begin(), rightHits.end(),
              [](CoverageCalculator& c1, CoverageCalculator& c2) -> bool {
                return c1.targetID < c2.targetID;
               });
    // Take the intersection of these two hit lists
    // Adopted from : http://en.cppreference.com/w/cpp/algorithm/set_intersection
    {
        auto leftIt = leftHits.begin();
        auto leftEnd = leftHits.end();
        auto rightIt = rightHits.begin();
        auto rightEnd = rightHits.end();
        while (leftIt != leftEnd && rightIt != rightEnd) {
            if (leftIt->targetID < rightIt->targetID) {
                ++leftIt;
            } else {
                if (!(rightIt->targetID < leftIt->targetID)) {
                    jointHits.push_back({leftIt->targetID,
                                         static_cast<size_t>(std::distance(leftHits.begin(), leftIt)),
                                         static_cast<size_t>(std::distance(rightHits.begin(), rightIt))});
                    ++leftIt;
                }
                ++rightIt;
            }
        }
    }
    // End vector-based code

    /* map based code
    {
        auto notFound = maxHitList.end();
        for (auto& kv : minHitList) {
            uint64_t refID = kv.first;
            if (maxHitList.find(refID) != notFound) {
                jointHits.emplace_back(refID);
            }
        }
    }
    */

    // Check if the fragment generated orphaned
    // lightweight alignments.
    bool isOrphan = (jointHits.size() == 0);

    uint32_t firstTranscriptID = std::numeric_limits<uint32_t>::max();
    double bestScore = -std::numeric_limits<double>::max();
    bool sortedByTranscript = true;
    int32_t lastTranscriptId = std::numeric_limits<int32_t>::min();

    if (BOOST_UNLIKELY(isOrphan and allowOrphans)) {
        //std::vector<CoverageCalculator> allHits;
        //allHits.reserve(totalHits);
        bool foundValidHit{false};

        // search for a hit on the left
        for (auto& tHitList : leftHits) {
            auto transcriptID = tHitList.targetID;
            auto& covChain = tHitList;
            Transcript& t = transcripts[transcriptID];
            if (!t.hasAnchorFragment()) { continue; }

            covChain.computeBestChain(t, frag.first.seq);
            double score = covChain.bestHitScore;

    	    // make sure orphaned fragment is near the end of the transcript
	    	// if (!nearEndOfTranscript(covChain, t, 1000)) { continue; }

            if (score >= fOpt * bestScore and score >= cutoffLeft) {
                foundValidHit = true;

        		if (score > bestScore) { bestScore = score; }
                bool isForward = covChain.isForward();
                int32_t hitPos = covChain.bestHitPos;
                auto fmt = salmon::utils::hitType(hitPos, isForward);

                if (leftHitCount == 0) {
                    firstTranscriptID = transcriptID;
                } else if (hitList.isUniquelyMapped() and transcriptID != firstTranscriptID) {
                    hitList.isUniquelyMapped() = false;
                }

                if (transcriptID  < lastTranscriptId) {
                    sortedByTranscript = false;
                }

                alnList.emplace_back(transcriptID, fmt, score, hitPos);
                readHits += score;
                ++hitListCount;
                ++leftHitCount;
            }
        }

        // search for a hit on the right
        for (auto& tHitList : rightHits) {
            // Prior
            // auto transcriptID = tHitList.first;
            auto transcriptID = tHitList.targetID;
            auto& covChain = tHitList;
            Transcript& t = transcripts[transcriptID];
            if (!t.hasAnchorFragment()) { continue; }

            covChain.computeBestChain(t, frag.second.seq);
            double score = covChain.bestHitScore;

            // make sure orphaned fragment is near the end of the transcript
            // if (!nearEndOfTranscript(covChain, t, 1000)) { continue; }

            if (score >= fOpt * bestScore and score >= cutoffRight) {
                if (score > bestScore) { bestScore = score; }
                foundValidHit = true;
                bool isForward = covChain.isForward();
                int32_t hitPos = covChain.bestHitPos;
                auto fmt = salmon::utils::hitType(hitPos, isForward);
                if (leftHitCount == 0) {
                    firstTranscriptID = transcriptID;
                } else if (hitList.isUniquelyMapped() and transcriptID != firstTranscriptID) {
                    hitList.isUniquelyMapped() = false;
                }

                alnList.emplace_back(transcriptID, fmt, score, hitPos);
                readHits += score;
                ++hitListCount;
                ++leftHitCount;
            }
        }

        if (alnList.size() > 0) {
            auto newEnd = std::stable_partition(alnList.begin(), alnList.end(),
                           [bestScore, fOpt](SMEMAlignment& aln) -> bool {
                                return aln.score() >= fOpt * bestScore;
                           });
            alnList.resize(std::distance(alnList.begin(), newEnd));
            if (!sortedByTranscript) {
                std::sort(alnList.begin(), alnList.end(),
                          [](const SMEMAlignment& x, const SMEMAlignment& y) -> bool {
                           return x.transcriptID() < y.transcriptID();
                          });
            }
        } else {
            return;
            // If we didn't have any *significant* hits --- add any *trivial* orphan hits
            size_t totalHits = leftHits.size() + rightHits.size();
            std::vector<uint32_t> txpIDs;
            txpIDs.reserve(totalHits);
            std::vector<double> auxProbs;
            auxProbs.reserve(totalHits);

            size_t txpIDsHash{0};
            std::vector<CoverageCalculator> allHits;
            allHits.reserve(totalHits);
            std::merge(leftHits.begin(), leftHits.end(),
                       rightHits.begin(), rightHits.end(),
                       std::back_inserter(allHits),
                       [](CoverageCalculator& c1, CoverageCalculator& c2) -> bool {
                        return c1.targetID < c2.targetID;
                       });
            double totProb{0.0};
            for (auto& h : allHits) {
                boost::hash_combine(txpIDsHash, h.targetID);
                txpIDs.push_back(h.targetID);
                double refLen =  std::max(1.0, static_cast<double>(transcripts[h.targetID].RefLength));
                double startProb = 1.0 / refLen;
                auxProbs.push_back(startProb);
                totProb += startProb;
            }
            if (totProb > 0.0) {
                double norm = 1.0 / totProb;
                for (auto& p : auxProbs) { p *= norm; }

                TranscriptGroup tg(txpIDs, txpIDsHash);
                eqBuilder.addGroup(std::move(tg), auxProbs);
            } else {
                salmonOpts.jointLog->warn("Unexpected empty hit group [orphaned]");
            }
        }
    } else { // Not an orphan
        for (auto jhp : jointHits) {
            auto& jointHitPtr = jhp;
            auto transcriptID = jhp.transcriptID;
            Transcript& t = transcripts[transcriptID];
            auto& leftHitList = leftHits[jhp.leftIndex];
            leftHitList.computeBestChain(t, frag.first.seq);
            if (leftHitList.bestHitScore >= cutoffLeft) {
                auto& rightHitList = rightHits[jhp.rightIndex];

                rightHitList.computeBestChain(t, frag.second.seq);
                if (rightHitList.bestHitScore < cutoffRight) { continue; }

                auto end1Start = leftHitList.bestHitPos;
                auto end2Start = rightHitList.bestHitPos;

                double score = (leftHitList.bestHitScore + rightHitList.bestHitScore) * 0.5;
                if (score < fOpt * bestScore) { continue; }

                if (score > bestScore) {
                    bestScore = score;
                }

                uint32_t fragLength = std::abs(static_cast<int32_t>(end1Start) -
                                               static_cast<int32_t>(end2Start)) + rightReadLength;

                bool end1IsForward = leftHitList.isForward();
                bool end2IsForward = rightHitList.isForward();

                uint32_t end1Pos = (end1IsForward) ? leftHitList.bestHitPos : leftHitList.bestHitPos + leftReadLength;
                uint32_t end2Pos = (end2IsForward) ? rightHitList.bestHitPos : rightHitList.bestHitPos + rightReadLength;
        		bool canDovetail = false;
                auto fmt = salmon::utils::hitType(end1Pos, end1IsForward, leftReadLength, end2Pos, end2IsForward, rightReadLength, canDovetail);

                if (readHits == 0) {
                    firstTranscriptID = transcriptID;
                } else if (hitList.isUniquelyMapped() and transcriptID != firstTranscriptID) {
                     hitList.isUniquelyMapped() = false;
                }

                int32_t minHitPos = std::min(end1Pos, end2Pos);
                if (transcriptID  < lastTranscriptId) {
                    sortedByTranscript = false;
                }
                // ANCHOR TEST
                t.setAnchorFragment();
                alnList.emplace_back(transcriptID, fmt, score, minHitPos, fragLength);
                ++readHits;
                ++hitListCount;
            }
        } // end for jointHits
        if (alnList.size() > 0) {
            auto newEnd = std::stable_partition(alnList.begin(), alnList.end(),
                           [bestScore, fOpt](SMEMAlignment& aln) -> bool {
                                return aln.score() >= fOpt * bestScore;
                           });
            alnList.resize(std::distance(alnList.begin(), newEnd));
            if (!sortedByTranscript) {
                std::sort(alnList.begin(), alnList.end(),
                          [](const SMEMAlignment& x, const SMEMAlignment& y) -> bool {
                           return x.transcriptID() < y.transcriptID();
                          });
            }
        } else {
            // If we didn't have any *significant* hits --- add any *trivial* joint hits
            return;
            std::vector<uint32_t> txpIDs;
            txpIDs.reserve(jointHits.size());
            std::vector<double> auxProbs;
            auxProbs.reserve(jointHits.size());

            size_t txpIDsHash{0};
            double totProb{0.0};
            for (auto& h : jointHits) {
                boost::hash_combine(txpIDsHash, h.transcriptID);
                txpIDs.push_back(h.transcriptID);
                double refLen =  std::max(1.0, static_cast<double>(transcripts[h.transcriptID].RefLength));
                double startProb = 1.0 / refLen;
                auxProbs.push_back(startProb);
                totProb += startProb;
            }
            if (totProb > 0.0) {
            double norm = 1.0 / totProb;
            for (auto& p : auxProbs) { p *= norm; }

            TranscriptGroup tg(txpIDs, txpIDsHash);
            eqBuilder.addGroup(std::move(tg), auxProbs);
            } else {
                salmonOpts.jointLog->warn("Unexpected empty hit group [paired]");
            }
        }

    } // end else
}

/**
  *   Get hits for single-end fragment
  *
  *
  */
template <typename CoverageCalculator>
inline void getHitsForFragment(jellyfish::header_sequence_qual& frag,
                        SalmonIndex* sidx,
                        smem_i *itr,
                        const bwtintv_v *a,
                        smem_aux_t* auxHits,
                        mem_opt_t* memOptions,
                        ReadExperiment& readExp,
                        const SalmonOpts& salmonOpts,
                        double coverageThresh,
                        uint64_t& upperBoundHits,
                        AlignmentGroup<SMEMAlignment>& hitList,
                        uint64_t& hitListCount,
                        std::vector<Transcript>& transcripts) {

    uint64_t leftHitCount{0};

    //std::unordered_map<uint64_t, CoverageCalculator> hits;
    std::vector<CoverageCalculator> hits;

    auto& eqBuilder = readExp.equivalenceClassBuilder();

    uint32_t readLength{0};

    //---------- get hits ----------------------//
    {
        std::string readStr   = frag.seq;
        uint32_t readLen      = frag.seq.size();

        readLength = readLen;

        for (int p = 0; p < readLen; ++p) {
            readStr[p] = nst_nt4_table[static_cast<int>(readStr[p])];
        }

        char* readPtr = const_cast<char*>(readStr.c_str());

        collectHitsForRead(sidx, a, auxHits,
                            memOptions,
                            salmonOpts,
                            reinterpret_cast<const uint8_t*>(readStr.c_str()),
                            readLen,
                            hits);

    }

    upperBoundHits += (hits.size() > 0) ? 1 : 0;

    int32_t lastTranscriptId = std::numeric_limits<int32_t>::min();
    bool sortedByTranscript{true};
    double fOpt{0.95};
    double bestScore = -std::numeric_limits<double>::max();

    size_t readHits{0};
    auto& alnList = hitList.alignments();
    hitList.isUniquelyMapped() = true;
    alnList.clear();

    uint32_t firstTranscriptID = std::numeric_limits<uint32_t>::max();
    double cutoff{ coverageThresh };//* readLength};
    for (auto& tHitList : hits) {
        // Prior
        // auto hitID = tHitList.first;
        // auto& covVec = tHitList.second;
        auto hitID = tHitList.targetID;
        auto& covVec = tHitList;

        // Coverage score
        Transcript& t = transcripts[hitID];
        covVec.computeBestChain(t, frag.seq);
        double score = covVec.bestHitScore;
        if (score >= fOpt * bestScore and covVec.bestHitScore >= cutoff) {

            bool isForward = covVec.isForward();
            if (score < fOpt * bestScore) { continue; }

        	if (score > bestScore) { bestScore = score; }

            auto hitPos = covVec.bestHitPos;
            auto fmt = salmon::utils::hitType(hitPos, isForward);

            if (leftHitCount == 0) {
                firstTranscriptID = hitID;
            } else if (hitList.isUniquelyMapped() and hitID != firstTranscriptID) {
                hitList.isUniquelyMapped() = false;
            }

            auto transcriptID = hitID;

            if (transcriptID  < lastTranscriptId) {
                sortedByTranscript = false;
            }

            alnList.emplace_back(transcriptID, fmt, score, hitPos);
            readHits += score;
            ++hitListCount;
            ++leftHitCount;
        }
    }
    if (alnList.size() > 0) {
        auto newEnd = std::stable_partition(alnList.begin(), alnList.end(),
                [bestScore, fOpt](SMEMAlignment& aln) -> bool {
                return aln.score() >= fOpt * bestScore;
                });
        alnList.resize(std::distance(alnList.begin(), newEnd));
        if (!sortedByTranscript) {
            std::sort(alnList.begin(), alnList.end(),
                    [](const SMEMAlignment& x, const SMEMAlignment& y) -> bool {
                     return x.transcriptID() < y.transcriptID();
                    });
        }
    }
    else {
        // If we didn't have any *significant* hits --- add any *trivial* joint hits
        return;
        std::vector<uint32_t> txpIDs;
        txpIDs.reserve(hits.size());
        double uniProb = 1.0 / hits.size();
        std::vector<double> auxProbs(hits.size(), uniProb);

        size_t txpIDsHash{0};
        for (auto& h : hits) {
            boost::hash_combine(txpIDsHash, h.targetID);
            txpIDs.push_back(h.targetID);
        }

        TranscriptGroup tg(txpIDs, txpIDsHash);
        eqBuilder.addGroup(std::move(tg), auxProbs);
    }


}

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename ParserT, typename CoverageCalculator>
void processReadsMEM(ParserT* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               AlnGroupVec<QuasiAlignment>& structureVec,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               std::atomic<uint64_t>& validHits,
               std::atomic<uint64_t>& upperBoundHits,
               SalmonIndex* sidx,
               std::vector<Transcript>& transcripts,
               ForgettingMassCalculator& fmCalc,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               mem_opt_t* memOptions,
               const SalmonOpts& salmonOpts,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn,
               volatile bool& writeToCache) {
    	// ERROR
	salmonOpts.jointLog->error("Quasimapping cannot be used with the FMD index --- please report this bug on GitHub");
	std::exit(1);
}

template <typename ParserT, typename CoverageCalculator>
void processReadsMEM(ParserT* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               AlnGroupVec<SMEMAlignment>& structureVec,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               std::atomic<uint64_t>& validHits,
               std::atomic<uint64_t>& upperBoundHits,
               SalmonIndex* sidx,
               std::vector<Transcript>& transcripts,
               ForgettingMassCalculator& fmCalc,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               mem_opt_t* memOptions,
               const SalmonOpts& salmonOpts,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn,
               volatile bool& writeToCache) {
  uint64_t count_fwd = 0, count_bwd = 0;
  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  // Super-MEM iterator
  smem_i *itr = smem_itr_init(sidx->bwaIndex()->bwt);
  const bwtintv_v *a = nullptr;
  smem_aux_t* auxHits = smem_aux_init();

  auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};

  while(true) {
    typename ParserT::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    rangeSize = j->nb_filled;
    if (rangeSize > structureVec.size()) {
        salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} --- this shouldn't happen.\n"
                                   "Please report this bug on GitHub", rangeSize, structureVec.size());
        std::exit(1);
    }

    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
        localUpperBoundHits = 0;

        auto& hitList = structureVec[i];
        getHitsForFragment<CoverageCalculator>(j->data[i], sidx, itr, a,
                                               auxHits,
                                               memOptions,
                                               readExp,
                                               salmonOpts,
                                               coverageThresh,
                                               localUpperBoundHits,
                                               hitList, hitListCount,
                                               transcripts);
        if (initialRound) {
            upperBoundHits += localUpperBoundHits;
        }

        // If the read mapped to > maxReadOccs places, discard it
        if (hitList.size() > salmonOpts.maxReadOccs ) { hitList.alignments().clear(); }
        validHits += hitList.size();
        locRead++;
        ++numObservedFragments;
        if (numObservedFragments % 50000 == 0) {
    	    iomutex.lock();
            const char RESET_COLOR[] = "\x1b[0m";
            char green[] = "\x1b[30m";
            green[3] = '0' + static_cast<char>(fmt::GREEN);
            char red[] = "\x1b[30m";
            red[3] = '0' + static_cast<char>(fmt::RED);
            if (initialRound) {
                fmt::print(stderr, "\033[A\r\r{}processed{} {} {}fragments{}\n", green, red, numObservedFragments, green, RESET_COLOR);
                fmt::print(stderr, "hits per frag:  {}; hit upper bound: {}",
                           validHits / static_cast<float>(prevObservedFrags),
                           upperBoundHits.load());
            } else {
                fmt::print(stderr, "\r\r{}processed{} {} {}fragments{}", green, red, numObservedFragments, green, RESET_COLOR);
            }
    	    iomutex.unlock();
        }


    } // end for i < j->nb_filled

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<SMEMAlignment> hitLists = boost::make_iterator_range(structureVec.begin(), structureVec.begin() + rangeSize);
    processMiniBatch<SMEMAlignment>(readExp, fmCalc,firstTimestepOfRound, rl, salmonOpts, hitLists, transcripts, clusterForest,
                     fragLengthDist, numAssignedFragments, eng, initialRound, burnedIn);
  }
  smem_aux_destroy(auxHits);
  smem_itr_destroy(itr);
}



/// START QUASI

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename RapMapIndexT>
void processReadsQuasi(paired_parser* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               AlnGroupVec<SMEMAlignment>& structureVec,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               std::atomic<uint64_t>& validHits,
               std::atomic<uint64_t>& upperBoundHits,
               RapMapIndexT* idx,
               std::vector<Transcript>& transcripts,
               ForgettingMassCalculator& fmCalc,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               mem_opt_t* memOptions,
               const SalmonOpts& salmonOpts,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn,
               volatile bool& writeToCache) {

    	// ERROR
	salmonOpts.jointLog->error("MEM-mapping cannot be used with the Quasi index --- please report this bug on GitHub");
	std::exit(1);
}

template <typename RapMapIndexT>
void processReadsQuasi(single_parser* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               AlnGroupVec<SMEMAlignment>& structureVec,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               std::atomic<uint64_t>& validHits,
               std::atomic<uint64_t>& upperBoundHits,
               RapMapIndexT* sidx,
               std::vector<Transcript>& transcripts,
               ForgettingMassCalculator& fmCalc,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               mem_opt_t* memOptions,
               const SalmonOpts& salmonOpts,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn,
               volatile bool& writeToCache) {
    	// ERROR
	salmonOpts.jointLog->error("MEM-mapping cannot be used with the Quasi index --- please report this bug on GitHub");
	std::exit(1);
}

template <typename RapMapIndexT>
void processReadsQuasi(paired_parser* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               AlnGroupVec<QuasiAlignment>& structureVec,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               std::atomic<uint64_t>& validHits,
               std::atomic<uint64_t>& upperBoundHits,
               RapMapIndexT* qidx,
               std::vector<Transcript>& transcripts,
               ForgettingMassCalculator& fmCalc,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               mem_opt_t* memOptions,
               const SalmonOpts& salmonOpts,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn,
               volatile bool& writeToCache) {
  uint64_t count_fwd = 0, count_bwd = 0;
  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  auto expectedLibType = rl.format();

  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};

  bool tooManyHits{false};
  size_t maxNumHits{salmonOpts.maxReadOccs};
  size_t readLen{0};
  SACollector<RapMapIndexT> hitCollector(qidx);
  SASearcher<RapMapIndexT> saSearcher(qidx);
  std::vector<QuasiAlignment> leftHits;
  std::vector<QuasiAlignment> rightHits;
  rapmap::utils::HitCounters hctr;

  while(true) {
    typename paired_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    rangeSize = j->nb_filled;
    if (rangeSize > structureVec.size()) {
        salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} --- this shouldn't happen.\n"
                                   "Please report this bug on GitHub", rangeSize, structureVec.size());
        std::exit(1);
    }

    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
        readLen = j->data[i].first.seq.length();
        tooManyHits = false;
        localUpperBoundHits = 0;
        auto& jointHitGroup = structureVec[i];
        auto& jointHits = jointHitGroup.alignments();
        jointHitGroup.clearAlignments();
        leftHits.clear();
        rightHits.clear();

        bool lh = hitCollector(j->data[i].first.seq,
                               leftHits, saSearcher,
                               MateStatus::PAIRED_END_LEFT);
        bool rh = hitCollector(j->data[i].second.seq,
                               rightHits, saSearcher,
                               MateStatus::PAIRED_END_RIGHT);

        rapmap::utils::mergeLeftRightHits(
                               leftHits, rightHits, jointHits,
                               readLen, maxNumHits, tooManyHits, hctr);

        if (initialRound) {
            upperBoundHits += (jointHits.size() > 0);
        }

        // If the read mapped to > maxReadOccs places, discard it
        if (jointHits.size() > salmonOpts.maxReadOccs ) { jointHitGroup.clearAlignments(); }


		// If we have mappings, then process them.
		if (jointHits.size() > 0) {
			bool isPaired = jointHits.front().mateStatus == rapmap::utils::MateStatus::PAIRED_END_PAIRED;
			// If we are ignoring orphans
			if (!salmonOpts.allowOrphans) {
				// If the mappings for the current read are not properly-paired (i.e. are orphans)
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
							return q.mateStatus == rapmap::utils::MateStatus::PAIRED_END_LEFT;
							});
					// Merge the hits so that the entire list is in order
					// by transcript ID.
					std::inplace_merge(jointHits.begin(), leftHitEndIt, jointHits.end(),
							[](const QuasiAlignment& a, const QuasiAlignment& b) -> bool {
							return a.transcriptID() < b.transcriptID();
							});
				}
			}

			for (auto& h : jointHits) {
				switch (h.mateStatus) {
					case MateStatus::PAIRED_END_LEFT:
						{
							h.format = salmon::utils::hitType(h.pos, h.fwd);
						}
						break;
					case MateStatus::PAIRED_END_RIGHT:
						{
							h.format = salmon::utils::hitType(h.pos, h.fwd);
						}
						break;
					case MateStatus::PAIRED_END_PAIRED:
						{
							uint32_t end1Pos = (h.fwd) ? h.pos : h.pos + h.readLen;
							uint32_t end2Pos = (h.mateIsFwd) ? h.matePos : h.matePos + h.mateLen;
							bool canDovetail = false;
							h.format = salmon::utils::hitType(end1Pos, h.fwd, h.readLen,
									end2Pos, h.mateIsFwd, h.mateLen, canDovetail);
						}
						break;
				}
			}
		} // If we have no mappings --- then there's nothing to do

        validHits += jointHits.size();
        locRead++;
        ++numObservedFragments;
        if (numObservedFragments % 500000 == 0) {
    	    iomutex.lock();
            const char RESET_COLOR[] = "\x1b[0m";
            char green[] = "\x1b[30m";
            green[3] = '0' + static_cast<char>(fmt::GREEN);
            char red[] = "\x1b[30m";
            red[3] = '0' + static_cast<char>(fmt::RED);
            if (initialRound) {
                fmt::print(stderr, "\033[A\r\r{}processed{} {} {}fragments{}\n", green, red, numObservedFragments, green, RESET_COLOR);
                fmt::print(stderr, "hits per frag:  {}; hit upper bound: {}",
                           validHits / static_cast<float>(prevObservedFrags),
                           upperBoundHits.load());
            } else {
                fmt::print(stderr, "\r\r{}processed{} {} {}fragments{}", green, red, numObservedFragments, green, RESET_COLOR);
            }
    	    iomutex.unlock();
        }


    } // end for i < j->nb_filled

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<QuasiAlignment> hitLists = boost::make_iterator_range(structureVec.begin(), structureVec.begin() + rangeSize);
    processMiniBatch<QuasiAlignment>(readExp, fmCalc,firstTimestepOfRound, rl, salmonOpts, hitLists, transcripts, clusterForest,
                     fragLengthDist, numAssignedFragments, eng, initialRound, burnedIn);
  }
}

// SINGLE END

// To use the parser in the following, we get "jobs" until none is
// available. A job behaves like a pointer to the type
// jellyfish::sequence_list (see whole_sequence_parser.hpp).
template <typename RapMapIndexT>
void processReadsQuasi(single_parser* parser,
               ReadExperiment& readExp,
               ReadLibrary& rl,
               AlnGroupVec<QuasiAlignment>& structureVec,
               std::atomic<uint64_t>& numObservedFragments,
               std::atomic<uint64_t>& numAssignedFragments,
               std::atomic<uint64_t>& validHits,
               std::atomic<uint64_t>& upperBoundHits,
               RapMapIndexT* qidx,
               std::vector<Transcript>& transcripts,
               ForgettingMassCalculator& fmCalc,
               ClusterForest& clusterForest,
               FragmentLengthDistribution& fragLengthDist,
               mem_opt_t* memOptions,
               const SalmonOpts& salmonOpts,
               double coverageThresh,
	           std::mutex& iomutex,
               bool initialRound,
               bool& burnedIn,
               volatile bool& writeToCache) {
  uint64_t count_fwd = 0, count_bwd = 0;
  // Seed with a real random value, if available
  std::random_device rd;

  // Create a random uniform distribution
  std::default_random_engine eng(rd());

  uint64_t prevObservedFrags{1};
  uint64_t leftHitCount{0};
  uint64_t hitListCount{0};

  auto expectedLibType = rl.format();


  uint64_t firstTimestepOfRound = fmCalc.getCurrentTimestep();

  size_t locRead{0};
  uint64_t localUpperBoundHits{0};
  size_t rangeSize{0};

  bool tooManyHits{false};
  size_t readLen{0};
  size_t maxNumHits{salmonOpts.maxReadOccs};
  SACollector<RapMapIndexT> hitCollector(qidx);
  SASearcher<RapMapIndexT> saSearcher(qidx);
  rapmap::utils::HitCounters hctr;

  while(true) {
    typename single_parser::job j(*parser); // Get a job from the parser: a bunch of read (at most max_read_group)
    if(j.is_empty()) break;           // If got nothing, quit

    rangeSize = j->nb_filled;
    if (rangeSize > structureVec.size()) {
        salmonOpts.jointLog->error("rangeSize = {}, but structureVec.size() = {} --- this shouldn't happen.\n"
                                   "Please report this bug on GitHub", rangeSize, structureVec.size());
        std::exit(1);
    }

    for(size_t i = 0; i < j->nb_filled; ++i) { // For all the read in this batch
        readLen = j->data[i].seq.length();
        tooManyHits = false;
        localUpperBoundHits = 0;
        auto& jointHitGroup = structureVec[i];
        auto& jointHits = jointHitGroup.alignments();
        jointHitGroup.clearAlignments();

        bool lh = hitCollector(j->data[i].seq,
                               jointHits, saSearcher,
                               MateStatus::SINGLE_END);

        if (initialRound) {
            upperBoundHits += (jointHits.size() > 0);
        }

        // If the read mapped to > maxReadOccs places, discard it
        if (jointHits.size() > salmonOpts.maxReadOccs ) { jointHitGroup.clearAlignments(); }

        for (auto& h : jointHits) {
            switch (h.mateStatus) {
                case MateStatus::SINGLE_END:
                    {
                        h.format = salmon::utils::hitType(h.pos, h.fwd);
                    }
                    break;
            }
        }

        validHits += jointHits.size();
        locRead++;
        ++numObservedFragments;
        if (numObservedFragments % 500000 == 0) {
    	    iomutex.lock();
            const char RESET_COLOR[] = "\x1b[0m";
            char green[] = "\x1b[30m";
            green[3] = '0' + static_cast<char>(fmt::GREEN);
            char red[] = "\x1b[30m";
            red[3] = '0' + static_cast<char>(fmt::RED);
            if (initialRound) {
                fmt::print(stderr, "\033[A\r\r{}processed{} {} {}fragments{}\n", green, red, numObservedFragments, green, RESET_COLOR);
                fmt::print(stderr, "hits per frag:  {}; hit upper bound: {}",
                           validHits / static_cast<float>(prevObservedFrags),
                           upperBoundHits.load());
            } else {
                fmt::print(stderr, "\r\r{}processed{} {} {}fragments{}", green, red, numObservedFragments, green, RESET_COLOR);
            }
    	    iomutex.unlock();
        }


    } // end for i < j->nb_filled

    prevObservedFrags = numObservedFragments;
    AlnGroupVecRange<QuasiAlignment> hitLists = boost::make_iterator_range(structureVec.begin(), structureVec.begin() + rangeSize);
    processMiniBatch<QuasiAlignment>(readExp, fmCalc,firstTimestepOfRound, rl, salmonOpts, hitLists, transcripts, clusterForest,
                     fragLengthDist, numAssignedFragments, eng, initialRound, burnedIn);
  }
}

/// DONE QUASI





int performBiasCorrection(boost::filesystem::path featPath,
                          boost::filesystem::path expPath,
                          double estimatedReadLength,
                          double kmersPerRead,
                          uint64_t mappedKmers,
                          uint32_t merLen,
                          boost::filesystem::path outPath,
                          size_t numThreads);

template <typename AlnT>
void processReadLibrary(
        ReadExperiment& readExp,
        ReadLibrary& rl,
        SalmonIndex* sidx,
        std::vector<Transcript>& transcripts,
        ClusterForest& clusterForest,
        std::atomic<uint64_t>& numObservedFragments, // total number of reads we've looked at
        std::atomic<uint64_t>& numAssignedFragments, // total number of assigned reads
        std::atomic<uint64_t>& upperBoundHits, // upper bound on # of mapped frags
        bool initialRound,
        bool& burnedIn,
        ForgettingMassCalculator& fmCalc,
        FragmentLengthDistribution& fragLengthDist,
        mem_opt_t* memOptions,
        const SalmonOpts& salmonOpts,
        double coverageThresh,
        bool greedyChain,
        std::mutex& iomutex,
        size_t numThreads,
        std::vector<AlnGroupVec<AlnT>>& structureVec,
        volatile bool& writeToCache){

            std::vector<std::thread> threads;

            std::atomic<uint64_t> numValidHits{0};
            rl.checkValid();

            auto indexType = sidx->indexType();

            std::unique_ptr<paired_parser> pairedParserPtr{nullptr};
            std::unique_ptr<single_parser> singleParserPtr{nullptr};
            // If the read library is paired-end
            // ------ Paired-end --------
            if (rl.format().type == ReadType::PAIRED_END) {


		    if (rl.mates1().size() != rl.mates2().size()) {
			    salmonOpts.jointLog->error("The number of provided files for "
					    "-1 and -2 must be the same!");
			    std::exit(1);
		    }

		    size_t numFiles = rl.mates1().size() + rl.mates2().size();
		    char** pairFileList = new char*[numFiles];
		    for (size_t i = 0; i < rl.mates1().size(); ++i) {
			    pairFileList[2*i] = const_cast<char*>(rl.mates1()[i].c_str());
			    pairFileList[2*i+1] = const_cast<char*>(rl.mates2()[i].c_str());
		    }

		    size_t maxReadGroup{miniBatchSize}; // Number of reads in each "job"
		    size_t concurrentFile{2}; // Number of files to read simultaneously
		    pairedParserPtr.reset(new
				    paired_parser(4 * numThreads, maxReadGroup,
					    concurrentFile, pairFileList, pairFileList+numFiles));

		    switch (indexType) {
			case SalmonIndexType::FMD:
			    {
				for(int i = 0; i < numThreads; ++i)  {
				    // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
				    // change value before the lambda below is evaluated --- crazy!
				    auto threadFun = [&,i]() -> void {
					processReadsMEM<paired_parser, TranscriptHitList>(
						pairedParserPtr.get(),
						readExp,
						rl,
						structureVec[i],
						numObservedFragments,
						numAssignedFragments,
						numValidHits,
						upperBoundHits,
						sidx,
						transcripts,
						fmCalc,
						clusterForest,
						fragLengthDist,
						memOptions,
						salmonOpts,
						coverageThresh,
						iomutex,
						initialRound,
						burnedIn,
						writeToCache);
				    };
				    threads.emplace_back(threadFun);
				}
				break;
				case SalmonIndexType::QUASI:
				{
            // True if we have a 64-bit SA index, false otherwise
            bool largeIndex = sidx->is64BitQuasi();
				    for(int i = 0; i < numThreads; ++i)  {
					// NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
					// change value before the lambda below is evaluated --- crazy!
          if (largeIndex) {
            auto threadFun = [&,i]() -> void {
              processReadsQuasi<RapMapSAIndex<int64_t>>(
                pairedParserPtr.get(),
                readExp,
                rl,
                structureVec[i],
                numObservedFragments,
                numAssignedFragments,
                numValidHits,
                upperBoundHits,
                sidx->quasiIndex64(),
                transcripts,
                fmCalc,
                clusterForest,
                fragLengthDist,
                memOptions,
                salmonOpts,
                coverageThresh,
                iomutex,
                initialRound,
                burnedIn,
                writeToCache);
              };
              threads.emplace_back(threadFun);
            } else {
              auto threadFun = [&,i]() -> void {
              processReadsQuasi<RapMapSAIndex<int32_t>>(
                pairedParserPtr.get(),
                readExp,
                rl,
                structureVec[i],
                numObservedFragments,
                numAssignedFragments,
                numValidHits,
                upperBoundHits,
                sidx->quasiIndex32(),
                transcripts,
                fmCalc,
                clusterForest,
                fragLengthDist,
                memOptions,
                salmonOpts,
                coverageThresh,
                iomutex,
                initialRound,
                burnedIn,
                writeToCache);
              };
              threads.emplace_back(threadFun);
            }

				    }
				}
				break;
			    } // end switch
		    }
		    for(int i = 0; i < numThreads; ++i) { threads[i].join(); }

            } // ------ Single-end --------
            else if (rl.format().type == ReadType::SINGLE_END) {

                char* readFiles[] = { const_cast<char*>(rl.unmated().front().c_str()) };
                size_t maxReadGroup{miniBatchSize}; // Number of files to read simultaneously
                size_t concurrentFile{1}; // Number of reads in each "job"
                stream_manager streams( rl.unmated().begin(),
                        rl.unmated().end(), concurrentFile);

                singleParserPtr.reset(new single_parser(4 * numThreads,
                                      maxReadGroup,
                                      concurrentFile,
                                      streams));

                switch (indexType) {
                    case SalmonIndexType::FMD:
                        {
                            for(int i = 0; i < numThreads; ++i)  {
                                // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
                                // change value before the lambda below is evaluated --- crazy!
                                auto threadFun = [&,i]() -> void {
                                    processReadsMEM<single_parser, TranscriptHitList>(
                                            singleParserPtr.get(),
                                            readExp,
                                            rl,
                                            structureVec[i],
                                            numObservedFragments,
                                            numAssignedFragments,
                                            numValidHits,
                                            upperBoundHits,
                                            sidx,
                                            transcripts,
                                            fmCalc,
                                            clusterForest,
                                            fragLengthDist,
                                            memOptions,
                                            salmonOpts,
                                            coverageThresh,
                                            iomutex,
                                            initialRound,
                                            burnedIn,
                                            writeToCache);
                                };
                                threads.emplace_back(threadFun);
                            }
                        }
                        break;

                    case SalmonIndexType::QUASI:
                    {
                      // True if we have a 64-bit SA index, false otherwise
                      bool largeIndex = sidx->is64BitQuasi();
                      for(int i = 0; i < numThreads; ++i)  {
                        // NOTE: we *must* capture i by value here, b/c it can (sometimes, does)
                        // change value before the lambda below is evaluated --- crazy!
                        if (largeIndex) {
                          auto threadFun = [&,i]() -> void {
                            processReadsQuasi<RapMapSAIndex<int64_t>>(
                              singleParserPtr.get(),
                              readExp,
                              rl,
                              structureVec[i],
                              numObservedFragments,
                              numAssignedFragments,
                              numValidHits,
                              upperBoundHits,
                              sidx->quasiIndex64(),
                              transcripts,
                              fmCalc,
                              clusterForest,
                              fragLengthDist,
                              memOptions,
                              salmonOpts,
                              coverageThresh,
                              iomutex,
                              initialRound,
                              burnedIn,
                              writeToCache);
                            };
                            threads.emplace_back(threadFun);
                          } else {
                            auto threadFun = [&,i]() -> void {
                              processReadsQuasi<RapMapSAIndex<int32_t>>(
                                singleParserPtr.get(),
                                readExp,
                                rl,
                                structureVec[i],
                                numObservedFragments,
                                numAssignedFragments,
                                numValidHits,
                                upperBoundHits,
                                sidx->quasiIndex32(),
                                transcripts,
                                fmCalc,
                                clusterForest,
                                fragLengthDist,
                                memOptions,
                                salmonOpts,
                                coverageThresh,
                                iomutex,
                                initialRound,
                                burnedIn,
                                writeToCache);
                              };
                              threads.emplace_back(threadFun);
                            }

                        }
                      }
                      break;
                }
                for(int i = 0; i < numThreads; ++i) { threads[i].join(); }
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
void quantifyLibrary(
        ReadExperiment& experiment,
        bool greedyChain,
        mem_opt_t* memOptions,
        SalmonOpts& salmonOpts,
        double coverageThresh,
        size_t numRequiredFragments,
        uint32_t numQuantThreads) {

    bool burnedIn{false};
    std::atomic<uint64_t> upperBoundHits{0};
    //ErrorModel errMod(1.00);
    auto& refs = experiment.transcripts();
    size_t numTranscripts = refs.size();
    // The *total* number of fragments observed so far (over all passes through the data).
    std::atomic<uint64_t> numObservedFragments{0};
    uint64_t prevNumObservedFragments{0};
    // The *total* number of fragments assigned so far (over all passes through the data).
    std::atomic<uint64_t> totalAssignedFragments{0};
    uint64_t prevNumAssignedFragments{0};

    auto jointLog = spdlog::get("jointLog");

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
            bool didReset = (salmonOpts.disableMappingCache) ?
                            (experiment.reset()) :
                            (experiment.softReset());

            if (!didReset) {
                std::string errmsg = fmt::sprintf(
                  "\n\n======== WARNING ========\n"
                  "One of the provided read files: [{}] "
                  "is not a regular file and therefore can't be read from "
                  "more than once.\n\n"
                  "We observed only {} mapping fragments when we wanted at least {}.\n\n"
                  "Please consider re-running Salmon with these reads "
                  "as a regular file!\n"
                  "NOTE: If you received this warning from salmon but did not "
                  "disable the mapping cache (--disableMappingCache), then there \n"
                  "was some other problem. Please make sure, e.g., that you have not "
                  "run out of disk space.\n"
                  "==========================\n\n",
                  experiment.readFilesAsString(), numObservedFragments, numRequiredFragments);
                jointLog->warn() << errmsg;
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
        auto processReadLibraryCallback =  [&](
                ReadLibrary& rl, SalmonIndex* sidx,
                std::vector<Transcript>& transcripts, ClusterForest& clusterForest,
                FragmentLengthDistribution& fragLengthDist,
                std::atomic<uint64_t>& numAssignedFragments,
                size_t numQuantThreads, bool& burnedIn) -> void  {

            processReadLibrary<AlnT>(experiment, rl, sidx, transcripts, clusterForest,
                    numObservedFragments, totalAssignedFragments, upperBoundHits,
                    initialRound, burnedIn, fmCalc, fragLengthDist,
                    memOptions, salmonOpts, coverageThresh, greedyChain,
                    ioMutex, numQuantThreads,
                    groupVec, writeToCache);

            numAssignedFragments = totalAssignedFragments - prevNumAssignedFragments;
            prevNumAssignedFragments = totalAssignedFragments;
        };

        // Process all of the reads
        experiment.processReads(numQuantThreads, processReadLibraryCallback);
        experiment.setNumObservedFragments(numObservedFragments);

        //EQCLASS
        bool done = experiment.equivalenceClassBuilder().finish();
        // skip the extra online rounds
        terminate = true;

        initialRound = false;
        ++roundNum;
        fmt::print(stderr, "\n# observed = {} / # required = {}\n",
                   numObservedFragments, numRequiredFragments);
        fmt::print(stderr, "hard # assigned = {} / # observed (this round) = {} : "
                           "upper bound assigned = {} \033[A\033[A",
                   experiment.numAssignedFragments(),
                   numObservedFragments - numPrevObservedFragments,
                   upperBoundHits);
        salmonOpts.fileLog->info("\nAt end of round {}\n"
                                   "==================\n"
                                   "Observed {} total fragments ({} in most recent round)\n",
                                   roundNum - 1,
                                   numObservedFragments,
                                   numObservedFragments - numPrevObservedFragments);
    }
    fmt::print(stderr, "\n\n\n\n");

    if (numObservedFragments <= prevNumObservedFragments) {
        jointLog->warn() << "Something seems to be wrong with the calculation "
            "of the mapping rate.  The recorded ratio is likely wrong.  Please "
            "file this as a bug report.\n";
    } else {
        double upperBoundMappingRate =
            upperBoundHits.load() /
            static_cast<double>(numObservedFragments.load());
        experiment.setNumObservedFragments(numObservedFragments - prevNumObservedFragments);
        experiment.setUpperBoundHits(upperBoundHits.load());
        experiment.setEffetiveMappingRate(upperBoundMappingRate);
    }

        jointLog->info("Mapping rate = {}\%\n",
                   experiment.effectiveMappingRate() * 100.0);
    jointLog->info("finished quantifyLibrary()");
}

int performBiasCorrectionSalmon(
        boost::filesystem::path featureFile,
        boost::filesystem::path expressionFile,
        boost::filesystem::path outputFile,
        size_t numThreads);

int salmonQuantify(int argc, char *argv[]) {
    using std::cerr;
    using std::vector;
    using std::string;
    namespace bfs = boost::filesystem;
    namespace po = boost::program_options;

    bool biasCorrect{false};
    bool optChain{false};
    size_t requiredObservations;

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
    generic.add_options()
    ("version,v", "print version string")
    ("help,h", "produce help message")
    ("index,i", po::value<string>()->required(), "Salmon index")
    ("libType,l", po::value<std::string>()->required(), "Format string describing the library type")
    ("unmatedReads,r", po::value<vector<string>>(&unmatedReadFiles)->multitoken(),
     "List of files containing unmated reads of (e.g. single-end reads)")
    ("mates1,1", po::value<vector<string>>(&mate1ReadFiles)->multitoken(),
        "File containing the #1 mates")
    ("mates2,2", po::value<vector<string>>(&mate2ReadFiles)->multitoken(),
        "File containing the #2 mates")
    ("allowOrphans", po::bool_switch(&(sopt.allowOrphans))->default_value(false), "Consider orphaned reads as valid hits when "
                        "performing lightweight-alignment.  This option will increase sensitivity (allow more reads to map and "
                        "more transcripts to be detected), but may decrease specificity as orphaned alignments are more likely "
                        "to be spurious -- this option is *always* set to true when using quasi-mapping.")
    ("threads,p", po::value<uint32_t>(&(sopt.numThreads))->default_value(sopt.numThreads), "The number of threads to use concurrently.")
    ("incompatPrior", po::value<double>(&(sopt.incompatPrior))->default_value(1e-20), "This option "
                        "sets the prior probability that an alignment that disagrees with the specified "
                        "library type (--libType) results from the true fragment origin.  Setting this to 0 "
                        "specifies that alignments that disagree with the library type should be \"impossible\", "
                        "while setting it to 1 says that alignments that disagree with the library type are no "
                        "less likely than those that do")
    ("minLen,k", po::value<int>(&(memOptions->min_seed_len))->default_value(19), "(S)MEMs smaller than this size won't be considered.")
    ("sensitive", po::bool_switch(&(sopt.sensitive))->default_value(false), "Setting this option enables the splitting of SMEMs that are larger "
                                        "than 1.5 times the minimum seed length (minLen/k above).  This may reveal high scoring chains of MEMs "
                                        "that are masked by long SMEMs.  However, this option makes lightweight-alignment a bit slower and is "
                                        "usually not necessary if the reference is of reasonable quality.")
    ("extraSensitive", po::bool_switch(&(sopt.extraSeedPass))->default_value(false), "Setting this option enables an extra pass of \"seed\" search. "
                                        "Enabling this option may improve sensitivity (the number of reads having sufficient coverage), but will "
                                        "typically slow down quantification by ~40%.  Consider enabling this option if you find the mapping rate to "
                                        "be significantly lower than expected.")
    ("coverage,c", po::value<double>(&coverageThresh)->default_value(0.70), "required coverage of read by union of SMEMs to consider it a \"hit\".")
    ("output,o", po::value<std::string>()->required(), "Output quantification file.")
    ("biasCorrect", po::value(&biasCorrect)->zero_tokens(), "[Experimental: Output both bias-corrected and non-bias-corrected "
                                                               "qunatification estimates.")
    ("geneMap,g", po::value<string>(), "File containing a mapping of transcripts to genes.  If this file is provided "
                                        "Salmon will output both quant.sf and quant.genes.sf files, where the latter "
                                        "contains aggregated gene-level abundance estimates.  The transcript to gene mapping "
                                        "should be provided as either a GTF file, or a in a simple tab-delimited format "
                                        "where each line contains the name of a transcript and the gene to which it belongs "
                                        "separated by a tab.  The extension of the file is used to determine how the file "
                                        "should be parsed.  Files ending in \'.gtf\' or \'.gff\' are assumed to be in GTF "
                                        "format; files with any other extension are assumed to be in the simple format.");
    //("optChain", po::bool_switch(&optChain)->default_value(false), "Chain MEMs optimally rather than greedily")

    // no sequence bias for now
    sopt.noSeqBiasModel = true;
    sopt.noRichEqClasses = false;
    // mapping cache has been deprecated
    sopt.disableMappingCache = true;

    po::options_description advanced("\n"
		    		     "advanced options");
    advanced.add_options()
    /*
    ("disableMappingCache", po::bool_switch(&(sopt.disableMappingCache))->default_value(false), "Setting this option disables the creation and use "
                                        "of the \"mapping cache\" file.  The mapping cache can speed up quantification significantly for smaller read "
                                        "libraries (i.e. where the number of mapped fragments is less than the required number of observations). However, "
                                        "for very large read libraries, the mapping cache is unnecessary, and disabling it may allow salmon to more effectively "
                                        "make use of a very large number of threads.")
    */
    ("fldMax" , po::value<size_t>(&(sopt.fragLenDistMax))->default_value(800), "The maximum fragment length to consider when building the empirical "
     											      "distribution")
    ("fldMean", po::value<size_t>(&(sopt.fragLenDistPriorMean))->default_value(200), "The mean used in the fragment length distribution prior")
    ("fldSD" , po::value<size_t>(&(sopt.fragLenDistPriorSD))->default_value(80), "The standard deviation used in the fragment length distribution prior")
    ("forgettingFactor,f", po::value<double>(&(sopt.forgettingFactor))->default_value(0.65), "The forgetting factor used "
                        "in the online learning schedule.  A smaller value results in quicker learning, but higher variance "
                        "and may be unstable.  A larger value results in slower learning but may be more stable.  Value should "
                        "be in the interval (0.5, 1.0].")
    ("maxOcc,m", po::value<int>(&(memOptions->max_occ))->default_value(200), "(S)MEMs occuring more than this many times won't be considered.")
    ("maxReadOcc,w", po::value<uint32_t>(&(sopt.maxReadOccs))->default_value(100), "Reads \"mapping\" to more than this many places won't be considered.")
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
    //("noSeqBiasModel", po::bool_switch(&(sopt.noSeqBiasModel))->default_value(false),
    //                    "Don't learn and apply a model of sequence-specific bias")
    ("numAuxModelSamples", po::value<uint32_t>(&(sopt.numBurninFrags))->default_value(5000000), "The first <numAuxModelSamples> are used to train the "
     			"auxiliary model parameters (e.g. fragment length distribution, bias, etc.).  After ther first <numAuxModelSamples> observations "
			"the auxiliary model parameters will be assumed to have converged and will be fixed.")
    ("numPreAuxModelSamples", po::value<uint32_t>(&(sopt.numPreBurninFrags))->default_value(1000000), "The first <numPreAuxModelSamples> will have their "
     			"assignment likelihoods and contributions to the transcript abundances computed without applying any auxiliary models.  The purpose "
			"of ignoring the auxiliary models for the first <numPreAuxModelSamples> observations is to avoid applying these models before thier "
			"parameters have been learned sufficiently well.")
    ("numRequiredObs,n", po::value(&requiredObservations)->default_value(50000000),
                                        "[Deprecated]: The minimum number of observations (mapped reads) that must be observed before "
                                        "the inference procedure will terminate.  If fewer mapped reads exist in the "
                                        "input file, then it will be read through multiple times.")
    ("splitWidth,s", po::value<int>(&(memOptions->split_width))->default_value(0), "If (S)MEM occurs fewer than this many times, search for smaller, contained MEMs. "
                                        "The default value will not split (S)MEMs, a higher value will result in more MEMs being explore and, thus, will "
                                        "result in increased running time.")
    ("splitSpanningSeeds,b", po::bool_switch(&(sopt.splitSpanningSeeds))->default_value(false), "Attempt to split seeds that happen to fall on the "
                                        "boundary between two transcripts.  This can improve the  fragment hit-rate, but is usually not necessary.")
    ("useVBOpt", po::bool_switch(&(sopt.useVBOpt))->default_value(false), "Use the Variational Bayesian EM rather than the "
     			"traditional EM algorithm for optimization in the batch passes.")
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
    all.add(generic).add(advanced).add(testing);

    po::options_description visible("salmon quant options");
    visible.add(generic).add(advanced);

    po::variables_map vm;
    try {
        auto orderedOptions = po::command_line_parser(argc,argv).
            options(all).run();

        po::store(orderedOptions, vm);

        if ( vm.count("help") ) {
            auto hstring = R"(
Quant
==========
Perform streaming mapping-based estimation of
transcript abundance from RNA-seq reads
)";
            std::cout << hstring << std::endl;
            std::cout << visible << std::endl;
            std::exit(1);
        }

        po::notify(vm);



        std::stringstream commentStream;
        commentStream << "# salmon (mapping-based) v" << salmon::version << "\n";
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

        // Verify the geneMap before we start doing any real work.
        bfs::path geneMapPath;
        if (vm.count("geneMap")) {
            // Make sure the provided file exists
            geneMapPath = vm["geneMap"].as<std::string>();
            if (!bfs::exists(geneMapPath)) {
                std::cerr << "Could not find transcript <=> gene map file " << geneMapPath << "\n";
                std::cerr << "Exiting now: please either omit the \'geneMap\' option or provide a valid file\n";
                std::exit(1);
            }
        }

        bool greedyChain = !optChain;
        bfs::path outputDirectory(vm["output"].as<std::string>());
        bfs::create_directory(outputDirectory);
        if (!(bfs::exists(outputDirectory) and bfs::is_directory(outputDirectory))) {
            std::cerr << "Couldn't create output directory " << outputDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }

        bfs::path indexDirectory(vm["index"].as<string>());
        bfs::path logDirectory = outputDirectory / "logs";

        sopt.indexDirectory = indexDirectory;
        sopt.outputDirectory = outputDirectory;

        // Create the logger and the logging directory
        bfs::create_directory(logDirectory);
        if (!(bfs::exists(logDirectory) and bfs::is_directory(logDirectory))) {
            std::cerr << "Couldn't create log directory " << logDirectory << "\n";
            std::cerr << "exiting\n";
            std::exit(1);
        }
        std::cerr << "Logs will be written to " << logDirectory.string() << "\n";

        bfs::path logPath = logDirectory / "salmon_quant.log";
	    // must be a power-of-two
        size_t max_q_size = 2097152;
        spdlog::set_async_mode(max_q_size);

        auto fileSink = std::make_shared<spdlog::sinks::simple_file_sink_mt>(logPath.string(), true);
        auto consoleSink = std::make_shared<spdlog::sinks::stderr_sink_mt>();
        auto consoleLog = spdlog::create("stderrLog", {consoleSink});
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

        {
            if (sopt.noFragLengthDist and !sopt.noEffectiveLengthCorrection) {
                jointLog->info() << "Error: You cannot enable --noFragLengthDist without "
                                 << "also enabling --noEffectiveLengthCorrection; exiting!\n";
                jointLog->flush();
                std::exit(1);
            }
        }

        // maybe arbitrary, but if it's smaller than this, consider it
        // equal to LOG_0
        if (sopt.incompatPrior < 1e-320) {
            sopt.incompatPrior = salmon::math::LOG_0;
        } else {
            sopt.incompatPrior = std::log(sopt.incompatPrior);
        }
        // END: option checking

        jointLog->info() << "parsing read library format";

        vector<ReadLibrary> readLibraries = salmon::utils::extractReadLibraries(orderedOptions);

        SalmonIndexVersionInfo versionInfo;
        boost::filesystem::path versionPath = indexDirectory / "versionInfo.json";
        versionInfo.load(versionPath);
        versionInfo.indexType();

        ReadExperiment experiment(readLibraries, indexDirectory, sopt);

        // Parameter validation
        // If we're allowing orphans, make sure that the read libraries are paired-end.
        // Otherwise, this option makes no sense.
        /*
        if (sopt.allowOrphans) {
            for (auto& rl : readLibraries) {
                if (!rl.isPairedEnd()) {
                    jointLog->error("You cannot specify the --allowOrphans argument "
                                    "for single-end libraries; exiting!");
                    std::exit(1);
                }
            }
        }
        */
        // end parameter validation


        // This will be the class in charge of maintaining our
    	// rich equivalence classes
        experiment.equivalenceClassBuilder().start();

        auto indexType = experiment.getIndex()->indexType();

        switch (indexType) {
            case SalmonIndexType::FMD:
                quantifyLibrary<SMEMAlignment>(experiment, greedyChain, memOptions, sopt, coverageThresh,
                                                requiredObservations, sopt.numThreads);
                break;
            case SalmonIndexType::QUASI:
                {
                    sopt.allowOrphans = true;
                    sopt.useQuasi = true;
                     quantifyLibrary<QuasiAlignment>(experiment, greedyChain, memOptions, sopt, coverageThresh,
                                                     requiredObservations, sopt.numThreads);
                }
                break;
        }

        // Now that the streaming pass is complete, we have
        // our initial estimates, and our rich equivalence
        // classes.  Perform further optimization until
        // convergence.
        CollapsedEMOptimizer optimizer;
        jointLog->info("Starting optimizer");
    	salmon::utils::normalizeAlphas(sopt, experiment);
        optimizer.optimize(experiment, sopt, 0.01, 10000);
        jointLog->info("Finished optimizer");

        free(memOptions);
        size_t tnum{0};

        jointLog->info("writing output \n");

        bfs::path estFilePath = outputDirectory / "quant.sf";

        commentStream << "# [ mapping rate ] => { " << experiment.effectiveMappingRate() * 100.0 << "\% }\n";
        commentString = commentStream.str();

        salmon::utils::writeAbundancesFromCollapsed(
                sopt, experiment, estFilePath, commentString);

        {
          bfs::path statPath = outputDirectory / "stats.tsv";
          std::ofstream statStream(statPath.string(), std::ofstream::out);
          statStream << "numObservedFragments\t" << experiment.numObservedFragsInFirstPass() << '\n';
          for (auto& t : experiment.transcripts()) {
              auto l = (sopt.noEffectiveLengthCorrection) ? t.RefLength : t.getCachedEffectiveLength();
              statStream << t.RefName << '\t' << l << '\n';
          }
          statStream.close();
        }

        if (sopt.numGibbsSamples > 0) {
            jointLog->info("Starting Gibbs Sampler");

            bfs::path gibbsSampleFile = sopt.outputDirectory / "quant_gibbs_samples.sf";
            sopt.jointLog->info("Writing posterior samples to {}", gibbsSampleFile.string());
            std::unique_ptr<BootstrapWriter> bsWriter(new TextBootstrapWriter(gibbsSampleFile, jointLog));
            bsWriter->writeHeader(commentString, experiment.transcripts());
            CollapsedGibbsSampler sampler;
            sampler.sample(experiment, sopt, bsWriter.get(), sopt.numGibbsSamples);

            jointLog->info("Finished Gibbs Sampler");
        } else if (sopt.numBootstraps > 0) {
            jointLog->info("Staring Bootstrapping");

            bfs::path bspath = outputDirectory / "quant_bootstraps.sf";
            std::unique_ptr<BootstrapWriter> bsWriter(new TextBootstrapWriter(bspath, jointLog));
            bsWriter->writeHeader(commentString, experiment.transcripts());
            optimizer.gatherBootstraps(experiment, sopt,
                    bsWriter.get(), 0.01, 10000);

            jointLog->info("Finished Bootstrapping");
        }


        // Now create a subdirectory for any parameters of interest
        bfs::path paramsDir = outputDirectory / "libParams";
        if (!boost::filesystem::exists(paramsDir)) {
            if (!boost::filesystem::create_directory(paramsDir)) {
                fmt::print(stderr, "{}ERROR{}: Could not create "
                           "output directory for experimental parameter "
                           "estimates [{}]. exiting.", ioutils::SET_RED,
                           ioutils::RESET_COLOR, paramsDir);
                std::exit(-1);
            }
        }

        bfs::path libCountFilePath = outputDirectory / "libFormatCounts.txt";
        experiment.summarizeLibraryTypeCounts(libCountFilePath);

        // Test writing out the fragment length distribution
        if (!sopt.noFragLengthDist) {
            bfs::path distFileName = paramsDir / "flenDist.txt";
            {
                std::unique_ptr<std::FILE, int (*)(std::FILE *)> distOut(std::fopen(distFileName.c_str(), "w"), std::fclose);
                fmt::print(distOut.get(), "{}\n", experiment.fragmentLengthDistribution()->toString());
            }
        }
        if (!sopt.noSeqBiasModel) {
            bfs::path biasFileName = paramsDir / "seqBias.txt";
            {
                std::unique_ptr<std::FILE, int (*)(std::FILE *)> biasOut(std::fopen(biasFileName.c_str(), "w"), std::fclose);
                fmt::print(biasOut.get(), "{}\n", experiment.sequenceBiasModel().toString());
            }
        }

        if (biasCorrect) {
            auto origExpressionFile = estFilePath;

            auto outputDirectory = estFilePath;
            outputDirectory.remove_filename();

            auto biasFeatPath = indexDirectory / "bias_feats.txt";
            auto biasCorrectedFile = outputDirectory / "quant_bias_corrected.sf";
            performBiasCorrectionSalmon(biasFeatPath, estFilePath, biasCorrectedFile, sopt.numThreads);
        }

        /** If the user requested gene-level abundances, then compute those now **/
        if (vm.count("geneMap")) {
            try {
                salmon::utils::generateGeneLevelEstimates(geneMapPath,
                                                            outputDirectory,
                                                            biasCorrect);
            } catch (std::invalid_argument& e) {
                fmt::print(stderr, "Error: [{}] when trying to compute gene-level "\
                                   "estimates. The gene-level file(s) may not exist",
                                   e.what());
            }
        }

    } catch (po::error &e) {
        std::cerr << "Exception : [" << e.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (const spdlog::spdlog_ex& ex) {
        std::cerr << "logger failed with : [" << ex.what() << "]. Exiting.\n";
        std::exit(1);
    } catch (std::exception& e) {
        std::cerr << "Exception : [" << e.what() << "]\n";
        std::cerr << argv[0] << " quant was invoked improperly.\n";
        std::cerr << "For usage information, try " << argv[0] << " quant --help\nExiting.\n";
        std::exit(1);
    }


    return 0;
}
