//I am not sure what are the things I would need soon
//let's go with what we have
#include <iostream>
#include <mutex>
#include <vector>
#include <random>
#include <unordered_map>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <thread>
#include <tuple>
#include <sstream>
#include <fstream>
#include <iostream>
#include <tuple>
#include <memory>
#include <cstring>
#include <queue>

//we already have timers

//cereal include
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <sparsepp/spp.h>
#include "parallel_hashmap/phmap.h"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/ostream_sink.h"
#include "spdlog/sinks/stdout_sinks.h"
#include "spdlog/sinks/ansicolor_sink.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/fmt/fmt.h"

#include "FastxParser.hpp"

//index header
#include "SelectiveAlignmentUtils.hpp"
#include "PuffAligner.hpp"
#include "ProgOpts.hpp"
#include "PufferfishIndex.hpp"
#include "PufferfishSparseIndex.hpp"
#include "PufferfishLossyIndex.hpp"
#include "Kmer.hpp"
#include "ScopedTimer.hpp"
#include "Util.hpp"
#include "SpinLock.hpp"
#include "MemCollector.hpp"
#include "SAMWriter.hpp"
#include "RefSeqConstructor.hpp"
#include "KSW2Aligner.hpp"
#include "zstr/zstr.hpp"


#define MATCH_SCORE 1
#define MISMATCH_SCORE -1
#define GAP_SCORE -1

#define EPS 5

#define ALLOW_VERBOSE 0

using paired_parser = fastx_parser::FastxParser<fastx_parser::ReadPair>;
using single_parser = fastx_parser::FastxParser<fastx_parser::ReadSeq>;

using HitCounters = pufferfish::util::HitCounters;
using QuasiAlignment = pufferfish::util::QuasiAlignment;
using MateStatus = pufferfish::util::MateStatus;

using MutexT = std::mutex;


//===========
// PAIRED END
//============
template<typename PufferfishIndexT>
void processReadsPair(paired_parser *parser,
                      PufferfishIndexT &pfi,
                      MutexT *iomutex,
                      std::shared_ptr<spdlog::logger> outQueue,
                      HitCounters &hctr,
                      phmap::flat_hash_set<std::string>& gene_names,
                      phmap::flat_hash_set<std::string>& rrna_names,
                      pufferfish::AlignmentOpts *mopts) {
    MemCollector<PufferfishIndexT> memCollector(&pfi);
    memCollector.configureMemClusterer(mopts->maxAllowedRefsPerHit);
    memCollector.setConsensusFraction(mopts->consensusFraction);
    memCollector.setAltSkip(mopts->altSkip);
    memCollector.setChainSubOptThresh(mopts->preMergeChainSubThresh);

    auto logger = spdlog::get("console");
    fmt::MemoryWriter sstream;
    BinWriter bstream;

    //size_t batchSize{2500} ;
    uint32_t readLen{0}, mateLen{0}, totLen{0};


    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> leftHits;
    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> rightHits;
    //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> leftHits;
    //phmap::flat_hash_map<size_t, std::vector<pufferfish::util::MemCluster>> rightHits;

    std::vector<pufferfish::util::MemCluster> recoveredHits;
    std::vector<pufferfish::util::JointMems> jointHits;
    PairedAlignmentFormatter<PufferfishIndexT *> formatter(&pfi);
    if (mopts->writeQualities) {
        formatter.enable_qualities();
    } else {
        formatter.disable_qualities();
    }
    pufferfish::util::QueryCache qc;

    //Initialize aligner ksw
    ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->mismatchScore);
    ksw2pp::KSW2Config config;

    config.dropoff = -1;
    config.gapo = mopts->gapOpenPenalty;
    config.gape = mopts->gapExtendPenalty;
    config.bandwidth = 15;
    config.flag = 0;
    config.flag |= KSW_EZ_RIGHT;
    config.flag |= KSW_EZ_SCORE_ONLY;
    aligner.config() = config;

    constexpr const int32_t invalidScore = std::numeric_limits<int32_t>::min();

    // don't think we should have these, they should
    // be baked into the index I think.
    bool filterGenomics = mopts->filterGenomics;
    bool filterMicrobiom = mopts->filterMicrobiom;
    bool filterBestScoreMicrobiom = mopts->filterMicrobiomBestScore;

    auto rg = parser->getReadGroup();

    phmap::flat_hash_map<uint32_t, std::pair<int32_t, int32_t>> bestScorePerTranscript;

    pufferfish::util::AlignmentConfig aconf;
    aconf.refExtendLength = mopts->refExtendLength;
    aconf.fullAlignment = mopts->fullAlignment;
    aconf.matchScore = mopts->matchScore;
    aconf.gapExtendPenalty = mopts->gapExtendPenalty;
    aconf.gapOpenPenalty = mopts->gapOpenPenalty;
    aconf.minScoreFraction = mopts->minScoreFraction;
    aconf.mimicBT2 = mopts->mimicBt2Default;
    aconf.allowOverhangSoftclip = mopts->allowOverhangSoftclip;
    aconf.allowSoftclip = mopts->allowSoftclip;
    aconf.alignmentMode = mopts->noOutput or !mopts->allowSoftclip ? pufferfish::util::PuffAlignmentMode::SCORE_ONLY : pufferfish::util::PuffAlignmentMode::APPROXIMATE_CIGAR;
    aconf.useAlignmentCache = mopts->useAlignmentCache;
    aconf.maxFragmentLength = mopts->maxFragmentLength;
    aconf.noDovetail = mopts->noDovetail;
    aconf.mismatchPenalty = mopts->mismatchScore;
    aconf.bestStrata = mopts->bestStrata;
    aconf.decoyPresent = mopts->filterGenomics or mopts->filterMicrobiom or mopts->filterMicrobiomBestScore;

    PuffAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(), aconf, aligner);

    std::vector<QuasiAlignment> jointAlignments;
    using pufferfish::util::BestHitReferenceType;
    BestHitReferenceType bestHitRefType{BestHitReferenceType::UNKNOWN};
    BestHitReferenceType hitRefType{BestHitReferenceType::UNKNOWN};

    pufferfish::util::MappingConstraintPolicy mpol;
    mpol.noDiscordant = mopts->noDiscordant;
    mpol.noOrphans = mopts->noOrphan;
    mpol.noDovetail = mopts->noDovetail;
    mpol.setPostMergeChainSubThresh(mopts->postMergeChainSubThresh);
    mpol.setOrphanChainSubThresh(mopts->orphanChainSubThresh);

    uint64_t firstDecoyIndex = pfi.firstDecoyIndex();

    //For filtering reads
    bool verbose = mopts->verbose;
//    auto &txpNames = pfi.getRefNames();
    uint32_t alignmentStreamLimit = mopts->alignmentStreamLimit;
    uint32_t alignmentStreamCount{0}, chunkReads{0};
    while (parser->refill(rg)) {
        for (auto read_it = rg.begin(); read_it != rg.end(); ++read_it) {
            auto& rpair = *read_it;
            readLen = static_cast<uint32_t >(rpair.first.seq.length());
            mateLen = static_cast<uint32_t >(rpair.second.seq.length());
            totLen = readLen + mateLen;
            bool tooShortRead = readLen < pfi.k();
            bool tooShortMate = mateLen < pfi.k();

            ++hctr.numReads;
            chunkReads++;
            jointHits.clear();
            leftHits.clear();
            rightHits.clear();
            recoveredHits.clear();
            memCollector.clear();
            jointAlignments.clear();

            // There is no way to revocer the following case other than aligning indels
            //verbose = rpair.first.seq == "CAGTGAGCCAAGATGGCGCCACTGCACTCCAGCCTGGGCAAAAAGAAACTCCATCTAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAAAACCCTGGTCCCT" or
            //          rpair.second.seq == "CAGTGAGCCAAGATGGCGCCACTGCACTCCAGCCTGGGCAAAAAGAAACTCCATCTAAAAAAAAAAAAAAAAAAAAAAAAAAGAGAAAACCCTGGTCCCT";

            // The only way to get it right is with non-heuristic chaining
            // verbose = rpair.first.seq == "AGCAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGTGGTGGGGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTAGAGAGGCACCAGCA" or
            //           rpair.second.seq == "AGCAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGTGGTGGGGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTGGTAGAGAGGCACCAGCA";

            //verbose = rpair.first.name == "mason_sample5_primary_1M_random.fasta.000050010/1";
            bool lh = tooShortRead ? false :
              memCollector(rpair.first.seq,
                           qc,
                           true, // isLeft
                           verbose);
            bool rh = tooShortMate ? false :
              memCollector(rpair.second.seq,
                           qc,
                           false, // isLeft
                           verbose);
            memCollector.findChains(rpair.first.seq,
                                   leftHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::PAIRED_END_LEFT,
                                   mopts->heuristicChaining,
                                   true, // isLeft
                                   verbose);
            memCollector.findChains(rpair.second.seq,
                                   rightHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::PAIRED_END_RIGHT,
                                   mopts->heuristicChaining,
                                   false, // isLeft
                                   verbose);

            hctr.numMappedAtLeastAKmer += (leftHits.size() > 0 || rightHits.size() > 0) ? 1 : 0;
            //do intersection on the basis of
            //performance, or going towards selective alignment
            //otherwise orphan

            // We also handle orphans inside this function
            /*if (rpair.first.name == "SRR3192367.14734") {
                std::stringstream ss;
                ss << "\n\nBefore joining:\n";
                for (auto &kv: leftHits) {
                    ss << kv.first << ":" << txpNames[kv.first] << " ";
                }
                ss << "\n";
                for (auto &kv: rightHits) {
                    ss << kv.first << ":" << txpNames[kv.first] << " ";
                }
                ss << "\n\n";
                std::cerr << ss.str();
            }*/
            if (tooShortRead and tooShortMate) {
              ++hctr.tooShortReads;
            } else {
              auto mergeRes = pufferfish::util::joinReadsAndFilter(leftHits, rightHits, jointHits,
                                                                 mopts->maxFragmentLength,
                                                                 totLen,
                                                                 mopts->scoreRatio,
                                                                 firstDecoyIndex,
                                                                 mpol, hctr);

              bool mergeStatusOR = (mergeRes == pufferfish::util::MergeResult::HAD_EMPTY_INTERSECTION or
                                  mergeRes == pufferfish::util::MergeResult::HAD_ONLY_LEFT or
                                  mergeRes == pufferfish::util::MergeResult::HAD_ONLY_RIGHT);

              if (mopts->recoverOrphans and mergeStatusOR and !tooShortRead and !tooShortMate) {
                // TODO NOTE : do futher testing
                bool recoveredAny = selective_alignment::utils::recoverOrphans(rpair.first.seq, rpair.second.seq, recoveredHits, jointHits, puffaligner, verbose);
                (void)recoveredAny;
              }

              hctr.peHits += jointHits.size();

#if ALLOW_VERBOSE
              if (verbose)
                 std::cerr<<"Number of hits: "<<jointHits.size()<<"\n";
#endif // ALLOW_VERBOSE
            }
            if (!mopts->justMap) {
              puffaligner.clear();
              int32_t bestScore = invalidScore;
              std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
              size_t idx{0};

                if (!mopts->genomicReads) { bestScorePerTranscript.clear(); }
                bestHitRefType = BestHitReferenceType::UNKNOWN;
                hitRefType = BestHitReferenceType::UNKNOWN;
                bool isMultimapping = (jointHits.size() > 1);
//                if (rpair.first.name == "SRR3192367.14734") {
//                    verbose = true;
//                }
//                std::stringstream ss;
//                if (verbose)
//                   ss << "\n\n found the read:\n" << rpair.first.name << " " << jointHits.size() <<"\n";
                puffaligner.getScoreStatus().reset();
                for (auto &&jointHit : jointHits) {
                    auto hitScore = puffaligner.calculateAlignments(rpair.first.seq, rpair.second.seq, jointHit, hctr, isMultimapping, verbose);
                    if (mopts->bestStrata and hitScore != invalidScore)
                        puffaligner.getScoreStatus().updateBest(hitScore - mopts->matchScore * std::max(rpair.first.seq.length(), rpair.second.seq.length()));
                    if ( (mopts->filterGenomics or mopts->filterMicrobiom or mopts->filterMicrobiomBestScore) and hitScore != invalidScore)
                        puffaligner.getScoreStatus().updateDecoy(hitScore - mopts->matchScore * std::max(rpair.first.seq.length(), rpair.second.seq.length()));
                    scores[idx] = hitScore;
//                    if (verbose)
//                        ss << txpNames[jointHit.tid] << " " << jointHit.alignmentScore << " " << scores[idx] << "\n";
                    const std::string& ref_name = pfi.refName(jointHit.tid);//txpNames[jointHit.tid];
                    if (filterMicrobiom and hitScore != invalidScore
                    and hitRefType != BestHitReferenceType::FILTERED) {
                        bool inFilteredList = (rrna_names.find(ref_name) != rrna_names.end());
                        if (inFilteredList) {
/*
                            std::stringstream ss;
                            ss << ref_name << " " << hitScore << "\n";
                            std::cerr << ss.str();
*/
                            hitRefType = BestHitReferenceType::FILTERED;
                        }
                    }
                    if (filterGenomics or filterBestScoreMicrobiom) {
                      bool inFilteredList = (gene_names.find(ref_name) != gene_names.end());
                      if (inFilteredList) {
//                          ss << "in filterMicrobiom\n";
                          scores[idx] = invalidScore;
                          if (hitScore > bestScore) {
                              bestHitRefType = BestHitReferenceType::FILTERED;
//                            if (verbose) ss << "1 " << bestScore;
                          } else if (hitScore == bestScore) {
                              bestHitRefType = (bestHitRefType == BestHitReferenceType::NON_FILTERED or
                                                bestHitRefType == BestHitReferenceType::BOTH) ?
                                               BestHitReferenceType::BOTH : BestHitReferenceType::FILTERED;
//                            if (verbose) ss << "2 " << bestScore;
                          }
                      } else { // the current hit was not in the list of references to filter
                          if (hitScore > bestScore) {
                              bestHitRefType = BestHitReferenceType::NON_FILTERED;
//                            if (verbose) ss << "3 " << bestScore;
                          } else if (hitScore == bestScore) {
                              bestHitRefType = (bestHitRefType == BestHitReferenceType::FILTERED or
                                                bestHitRefType == BestHitReferenceType::BOTH) ?
                                               BestHitReferenceType::BOTH : BestHitReferenceType::NON_FILTERED;
//                            if (verbose) ss << "4 " << bestScore;
                          }
                      }
//                      if (verbose) ss << "\n";
                    }

                    bestScore = (hitScore > bestScore) ? hitScore : bestScore;

                    // use ALIGNMENT SCORE, not COVERAGE
                    // We should have valid alignment scores at this point as we are inside !justMap if and passed calculating alignment score
                    if (!mopts->genomicReads) {
                        // removing dupplicate hits from a read to the same transcript
                      auto it = bestScorePerTranscript.find(jointHit.tid);
                        if (it == bestScorePerTranscript.end()) {
                          // if we didn't have any alignment for this transcript yet, then
                          // this is the current best
                          // cast is a hack :(
                          bestScorePerTranscript[jointHit.tid].first = jointHit.alignmentScore+jointHit.mateAlignmentScore;/*static_cast<int32_t>(jointHit.coverage());*/
                          bestScorePerTranscript[jointHit.tid].second = idx;
                        } else if (jointHit.coverage() > it->second.first) {
                          // otherwise, if we had an alignment for this transcript and it's
                          // better than the current best, then set the best score to this
                          // alignment's score, and invalidate the previous alignment
                          it->second.first = jointHit.alignmentScore+jointHit.mateAlignmentScore;/*static_cast<int32_t>(jointHit.coverage());*/
                          scores[it->second.second] = invalidScore;
                          it->second.second = idx;
                        } else {
                          // otherwise, there is already a better mapping for this transcript.
                          scores[idx] = invalidScore;
                        }
                    }
                    ++idx;
                }

                /*if (verbose) {
                    ss << "\n\nbestHitRefType: ";
                    switch(bestHitRefType) {
                        case BestHitReferenceType::UNKNOWN:
                            ss << "unknown\n"; break;
                        case BestHitReferenceType::FILTERED:
                            ss << "filtered\n"; break;
                        case BestHitReferenceType::NON_FILTERED:
                            ss << "nonFiltered\n"; break;
                        case BestHitReferenceType::BOTH:
                            ss << "both\n"; break;
                    }
                }
                if (verbose) {std::cerr << ss.str();  verbose = false; std::exit(1);}*/
                if (filterGenomics and (bestHitRefType == BestHitReferenceType::FILTERED) ) {
                    // This read is likely come from decoy sequence and should be discarded from reference alignments
                    continue;
                }
               /* std::stringstream ss;

                ss << filterMicrobiom << " ";
                switch(bestHitRefType) {
                    case BestHitReferenceType::UNKNOWN:
                        ss << "unknown\n"; break;
                    case BestHitReferenceType::FILTERED:
                        ss << "filtered\n"; break;
                    case BestHitReferenceType::NON_FILTERED:
                        ss << "nonFiltered\n"; break;
                    case BestHitReferenceType::BOTH:
                        ss << "both\n"; break;
                }
                std::cerr << ss.str();*/
                if (filterBestScoreMicrobiom and (bestHitRefType != BestHitReferenceType::NON_FILTERED)) {
                    continue;
                }
                if (filterMicrobiom and (hitRefType == BestHitReferenceType::FILTERED)) {
//                    std::cerr << "filtered\n";
                    continue;
                }
                // Filter out these alignments with low scores
                uint32_t ctr{0};
                if (bestScore > invalidScore) {
                  bool filterBestStrata = mopts->bestStrata;
                  jointHits.erase(
                                  std::remove_if(jointHits.begin(), jointHits.end(),
                                                 [&ctr, &scores, filterBestStrata, bestScore](pufferfish::util::JointMems &) -> bool {
                                                   bool rem = filterBestStrata ? (scores[ctr] < bestScore) : (scores[ctr] == invalidScore);
                                                   ++ctr;
                                                   return rem;
                                                 }),
                                  jointHits.end()
                                  );

                  if (mopts->primaryAlignment and !jointHits.empty()) {
                    jointHits.resize(1);
                  }
                } else {
                    // There is no alignment with high quality for this read, so we skip this reads' alignments
                    jointHits.clear();
                }
            }

            if (jointHits.size() > mopts->maxNumHits) {
                std::nth_element(jointHits.begin(), jointHits.begin() + mopts->maxNumHits,jointHits.end(),
                          [](const auto &lhs, const auto &rhs) {
                              return lhs.alignmentScore > rhs.alignmentScore;
                          });
                jointHits.erase(jointHits.begin() + mopts->maxNumHits, jointHits.end());
            }

            hctr.totHits += !jointHits.empty() && !jointHits.back().isOrphan() ? 1 : 0;;
            hctr.numMapped += !jointHits.empty() ? 1 : 0;
            if (mopts->noOrphan) {
                hctr.numOfOrphans += jointHits.empty() && (lh || rh);
            } else {
                hctr.numOfOrphans += !jointHits.empty() && (jointHits.back().isOrphan()) ? 1 : 0;
            }

            if (jointHits.size() > hctr.maxMultimapping) {
                hctr.maxMultimapping = jointHits.size();
            }

            for (auto &&jointHit : jointHits) {
              if (jointHit.isOrphan()) {
                readLen = jointHit.isLeftAvailable() ? readLen : mateLen;
                jointAlignments.emplace_back(jointHit.tid,           // reference id
                                             jointHit.orphanClust()->getTrFirstHitPos(),     // reference pos
                                             jointHit.orphanClust()->isFw,     // fwd direction
                                             readLen, // read length
                                             jointHit.orphanClust()->cigar, // cigar string
                                             jointHit.fragmentLen,       // fragment length
                                             false);
                auto &qaln = jointAlignments.back();
                // NOTE : score should not be filled in from a double
                qaln.score = mopts->justMap ? static_cast<int32_t >(jointHit.orphanClust()->coverage):jointHit.alignmentScore;//()->coverage;
                // NOTE : wth is numHits?
                qaln.numHits = static_cast<uint32_t >(jointHits.size());//orphanClust()->coverage;
                qaln.mateStatus = jointHit.mateStatus;
              } else {
                jointAlignments.emplace_back(jointHit.tid,           // reference id
                                             jointHit.leftClust->getTrFirstHitPos(),     // reference pos
                                             jointHit.leftClust->isFw,     // fwd direction
                                             readLen, // read length
                                             jointHit.leftClust->cigar, // cigar string
                                             jointHit.fragmentLen,       // fragment length
                                             true);         // properly paired
                // Fill in the mate info
                auto &qaln = jointAlignments.back();
                qaln.mateLen = mateLen;
                qaln.mateCigar = jointHit.rightClust->cigar;
                qaln.matePos = static_cast<int32_t >(jointHit.rightClust->getTrFirstHitPos());
                qaln.mateIsFwd = jointHit.rightClust->isFw;
                qaln.mateStatus = MateStatus::PAIRED_END_PAIRED;
                // NOTE : wth is numHits?
                qaln.numHits = static_cast<uint32_t >(jointHits.size());
                // NOTE : score should not be filled in from a double
                qaln.score = mopts->justMap ? static_cast<int32_t >(jointHit.leftClust->coverage):jointHit.alignmentScore;
                qaln.mateScore = mopts->justMap ? static_cast<int32_t >(jointHit.rightClust->coverage):jointHit.mateAlignmentScore;;
              }
            }

            hctr.totAlignment += jointAlignments.size();

            if (!mopts->noOutput) {
              if (mopts->krakOut) {
                writeAlignmentsToKrakenDump(rpair,  formatter,  jointHits, bstream, mopts->justMap);
                alignmentStreamCount += jointHits.size();
              } else if (mopts->salmonOut) {
                writeAlignmentsToKrakenDump(rpair,  formatter,  jointHits, bstream, mopts->justMap, false);
                alignmentStreamCount += jointHits.size();
              }  else if (mopts->radOut) {
                writeAlignmentsToRADPair(rpair,  formatter,  jointAlignments, bstream, false);
                alignmentStreamCount += jointAlignments.size();
              } else if (jointAlignments.size() > 0) {
                writeAlignmentsToStream(rpair, formatter, jointAlignments, sstream, !mopts->noOrphan);
                alignmentStreamCount += jointAlignments.size();
              } else if (jointAlignments.size() == 0) {
                writeUnalignedPairToStream(rpair, formatter, sstream);
                alignmentStreamCount += 1;
              }
            }

            // write them on cmd
            if (hctr.numReads > hctr.lastPrint + 100000) {
                hctr.lastPrint.store(hctr.numReads.load());
                if (!mopts->quiet and iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads : "
                              << "pe / read = " << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = " << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
                    iomutex->unlock();
                }
            }
            //puffaligner.clear();
            // try dumping the output
            bool last_read = (read_it + 1 == rg.end());
            if (!mopts->noOutput and (alignmentStreamCount > alignmentStreamLimit or last_read)) {
                // Get rid of last newline
                if (mopts->salmonOut) {
                    if (bstream.getBytes() != 0) {
                        BinWriter sbw(sizeof(uint64_t));
                        sbw << bstream.getBytes();
                        outQueue->info("{}{}", sbw, bstream);
                    }
                } else if (mopts->radOut) {
                    if (bstream.getBytes() != 0) {
                        BinWriter sbw(sizeof(uint64_t));
                        sbw << bstream.getBytes() << chunkReads;
                        outQueue->info("{}{}", sbw, bstream);
                        chunkReads = 0;
                    } 
                } else if (mopts->krakOut) {
                    outQueue->info("{}", bstream);
                } else {
                    std::string outStr(sstream.str());
                    if (!outStr.empty()) {
                        outStr.pop_back();
                        outQueue->info(std::move(outStr));
                    }
                }
                sstream.clear();
                bstream.clear();
                alignmentStreamCount = 0;
            }
        } // for all reads in this job
    } // processed all reads
}

//===========
// SINGLE END
//============
template<typename PufferfishIndexT>
void processReadsSingle(single_parser *parser,
                        PufferfishIndexT &pfi,
                        MutexT *iomutex,
                        std::shared_ptr<spdlog::logger> outQueue,
                        HitCounters &hctr,
                        phmap::flat_hash_set<std::string>& gene_names,
                        pufferfish::AlignmentOpts *mopts) {
    MemCollector<PufferfishIndexT> memCollector(&pfi);
    memCollector.configureMemClusterer(mopts->maxAllowedRefsPerHit);
    memCollector.setConsensusFraction(mopts->consensusFraction);
    memCollector.setAltSkip(mopts->altSkip);

    using pufferfish::util::BestHitReferenceType;
    BestHitReferenceType bestHitRefType{BestHitReferenceType::UNKNOWN};
    phmap::flat_hash_map<uint32_t, std::pair<int32_t, int32_t>> bestScorePerTranscript;

    auto logger = spdlog::get("console");
    fmt::MemoryWriter sstream;
    BinWriter bstream;
    //size_t batchSize{2500} ;
    uint32_t readLen{0};
    std::string dummyRead = "";
    //size_t totLen{0};

    pufferfish::util::CachedVectorMap<size_t, std::vector<pufferfish::util::MemCluster>, std::hash<size_t>> leftHits;
    std::vector<pufferfish::util::JointMems> jointHits;
    PairedAlignmentFormatter<PufferfishIndexT *> formatter(&pfi);
    if (mopts->writeQualities) {
        formatter.enable_qualities();
    } else {
        formatter.disable_qualities();
    }
    pufferfish::util::QueryCache qc;
    std::vector<pufferfish::util::MemCluster> all;

    //Initialize aligner ksw
    ksw2pp::KSW2Aligner aligner(mopts->matchScore, mopts->mismatchScore);
    ksw2pp::KSW2Config config;

    config.dropoff = -1;
    config.gapo = mopts->gapOpenPenalty;
    config.gape = mopts->gapExtendPenalty;
    config.bandwidth = 15;
    config.flag = 0;
    config.flag |= KSW_EZ_RIGHT;
    config.flag |= KSW_EZ_SCORE_ONLY;
    aligner.config() = config;

    constexpr const int32_t invalidScore = std::numeric_limits<int32_t>::min();

//    auto &txpNames = pfi.getRefNames();

    pufferfish::util::AlignmentConfig aconf;
    aconf.refExtendLength = mopts->refExtendLength;
    aconf.fullAlignment = mopts->fullAlignment;
    aconf.matchScore = mopts->matchScore;
    aconf.gapExtendPenalty = mopts->gapExtendPenalty;
    aconf.gapOpenPenalty = mopts->gapOpenPenalty;
    aconf.minScoreFraction = mopts->minScoreFraction;
    aconf.mimicBT2 = mopts->mimicBt2Default;
    aconf.allowOverhangSoftclip = mopts->allowOverhangSoftclip;
    aconf.allowSoftclip = mopts->allowSoftclip;
    aconf.alignmentMode = mopts->noOutput or !mopts->allowSoftclip ? pufferfish::util::PuffAlignmentMode::SCORE_ONLY : pufferfish::util::PuffAlignmentMode::APPROXIMATE_CIGAR;
    aconf.useAlignmentCache = mopts->useAlignmentCache;
    aconf.mismatchPenalty = mopts->mismatchScore;
    aconf.bestStrata = mopts->bestStrata;
    aconf.decoyPresent = mopts->filterGenomics or mopts->filterMicrobiom or mopts->filterMicrobiomBestScore;

    PuffAligner puffaligner(pfi.refseq_, pfi.refAccumLengths_, pfi.k(), aconf, aligner);

    uint32_t alignmentStreamLimit = mopts->alignmentStreamLimit;
    uint32_t alignmentStreamCount{0}, chunkReads{0};

    auto rg = parser->getReadGroup();
    while (parser->refill(rg)) {
        for (auto read_it = rg.begin(); read_it != rg.end(); ++read_it) {
            auto& read = *read_it;
            readLen = static_cast<uint32_t >(read.seq.length());
            bool tooShortRead = readLen < pfi.k();
            auto totLen = readLen;
            bool verbose = false;
            //if (verbose) std::cerr << read.name << "\n";
            ++hctr.numReads;

            jointHits.clear();
            leftHits.clear();
            memCollector.clear();


            bool filterGenomics = mopts->filterGenomics;
            bool filterMicrobiom = mopts->filterMicrobiom;

            bool lh = tooShortRead? false :
              memCollector(read.seq,
                           qc,
                           true, // isLeft
                           verbose);
            memCollector.findChains(read.seq,
                                   leftHits,
                                   mopts->maxSpliceGap,
                                   MateStatus::SINGLE_END,
                                   mopts->heuristicChaining,
                                   true, // isLeft
                                   verbose);

            (void) lh;
            all.clear();
            if (tooShortRead) {
              ++hctr.tooShortReads;
            } else {
              pufferfish::util::joinReadsAndFilterSingle(leftHits, jointHits,
                                     totLen,
                                     mopts->scoreRatio);
            }

            std::vector<QuasiAlignment> jointAlignments;
            std::vector<std::pair<uint32_t, std::vector<pufferfish::util::MemCluster>::iterator>> validHits;

            if (!mopts->justMap) {
                puffaligner.clear();

                int32_t bestScore = invalidScore;
                std::vector<decltype(bestScore)> scores(jointHits.size(), bestScore);
                size_t idx{0};

                bool bestScoreGenomic{false};
                bool bestScoreTxpomic{false};
                std::map<int32_t, std::vector<int32_t>> transcript_set;
                if (!mopts->genomicReads) { bestScorePerTranscript.clear(); }
                bestHitRefType = BestHitReferenceType::UNKNOWN;
                bool isMultimapping = (jointHits.size() > 1);
                puffaligner.getScoreStatus().reset();
                for (auto &jointHit : jointHits) {
                    int32_t hitScore = puffaligner.calculateAlignments(read.seq, jointHit, hctr, isMultimapping, verbose);
                    if (mopts->bestStrata) puffaligner.getScoreStatus().updateBest(hitScore);
                    if (mopts->filterGenomics or mopts->filterMicrobiom or mopts->filterMicrobiomBestScore) puffaligner.getScoreStatus().updateDecoy(hitScore);
                    scores[idx] = hitScore;

                    const std::string& ref_name = pfi.refName(jointHit.tid);//txpNames[jointHit.tid];
                    if (filterGenomics or filterMicrobiom) {
                        bool inFilteredList = (gene_names.find(ref_name) != gene_names.end());
                        if (inFilteredList) {
                            scores[idx] = invalidScore;
                            if (hitScore > bestScore) {
                                bestHitRefType = BestHitReferenceType::FILTERED;
                            } else if (hitScore == bestScore) {
                                bestHitRefType = (bestHitRefType == BestHitReferenceType::NON_FILTERED or
                                                  bestHitRefType == BestHitReferenceType::BOTH) ?
                                                 BestHitReferenceType::BOTH : BestHitReferenceType::FILTERED;
                            }
                        } else { // the current hit was not in the list of references to filter
                            if (hitScore > bestScore) {
                                bestHitRefType = BestHitReferenceType::NON_FILTERED;
                            } else if (hitScore == bestScore) {
                                bestHitRefType = (bestHitRefType == BestHitReferenceType::FILTERED or
                                                  bestHitRefType == BestHitReferenceType::BOTH) ?
                                                 BestHitReferenceType::BOTH : BestHitReferenceType::NON_FILTERED;
                            }
                        }
                    }

                    bestScore = (hitScore > bestScore) ? hitScore : bestScore;


                    // use ALIGNMENT SCORE, not COVERAGE
                    // We should have valid alignment scores at this point as we are inside !justMap if and passed calculating alignment score
                    if (!mopts->genomicReads) {
                        // removing dupplicate hits from a read to the same transcript
                        auto it = bestScorePerTranscript.find(jointHit.tid);
                        if (it == bestScorePerTranscript.end()) {
                            bestScorePerTranscript[jointHit.tid].first = jointHit.alignmentScore;
                            bestScorePerTranscript[jointHit.tid].second = idx;
                        } else if (jointHit.coverage() > it->second.first) {
                            it->second.first = jointHit.alignmentScore;
                            scores[it->second.second] = invalidScore;
                            it->second.second = idx;
                        } else {
                            // otherwise, there is already a better mapping for this transcript.
                            scores[idx] = invalidScore;
                        }
                    }
                    ++idx;
                }

                if (filterGenomics and bestScoreGenomic and !bestScoreTxpomic) {
                    // This read is likely come from the genome and should be discarded from txptomic alignments
                    continue;
                }
                if (filterMicrobiom and bestScoreGenomic) {
                    continue;
                }

                // Filter out these alignments with low scores
                uint32_t ctr{0};
                if (bestScore > invalidScore) {
                    bool filterBestStrata = mopts->bestStrata;
                    jointHits.erase(
                            std::remove_if(jointHits.begin(), jointHits.end(),
                                           [&ctr, &scores, filterBestStrata, bestScore](pufferfish::util::JointMems &) -> bool {
                                               bool rem = filterBestStrata ? (scores[ctr] < bestScore) : (scores[ctr] == invalidScore);
                                               ++ctr;
                                               return rem;
                                           }),
                            jointHits.end()
                    );

                    if (mopts->primaryAlignment and !jointHits.empty()) {
                        jointHits.resize(1);
                    }
                } else {
                    // There is no alignment with high quality for this read, so we skip this reads' alignments
                    jointHits.clear();
                }
            }

            if (jointHits.size() > mopts->maxNumHits) {
                std::nth_element(jointHits.begin(), jointHits.begin() + mopts->maxNumHits,jointHits.end(),
                          [](const auto &lhs, const auto &rhs) {
                              return lhs.alignmentScore > rhs.alignmentScore;
                          });
                jointHits.erase(jointHits.begin() + mopts->maxNumHits, jointHits.end());
            }

            hctr.numMappedAtLeastAKmer += jointHits.size() > 0 ? 1 : 0;
            hctr.totHits += jointHits.size();
            hctr.seHits += jointHits.size();
            hctr.numMapped += !jointHits.empty() ? 1 : 0;
            if (jointHits.size() > hctr.maxMultimapping) {
                hctr.maxMultimapping = jointHits.size();
            }

            for (auto &jointHit : jointHits) {
                jointAlignments.emplace_back(jointHit.tid,           // reference id
                                             jointHit.orphanClust()->getTrFirstHitPos(),     // reference pos
                                             jointHit.orphanClust()->isFw,     // fwd direction
                                             readLen, // read length
                                             jointHit.orphanClust()->cigar, // cigar string
                                             readLen,       // fragment length
                                             false);
                auto &qaln = jointAlignments.back();
                qaln.mateLen = readLen;
                qaln.mateCigar = "";
                qaln.matePos = 0;       // jointHit.rightClust->getTrFirstHitPos();
                qaln.mateIsFwd = false; // jointHit.rightClust->isFw;
                qaln.mateStatus = MateStatus::SINGLE_END;
                qaln.numHits = static_cast<uint32_t >(jointHits.size());//orphanClust()->coverage;
                qaln.score = jointHit.alignmentScore;
                qaln.mateScore = 0;
                if (mopts->justMap) jointHit.orphanClust()->coverage = jointHit.alignmentScore;
                validHits.emplace_back(jointHit.tid, jointHit.orphanClust());
            }

            hctr.totAlignment += jointHits.size();

            // write puffkrak format output
            if (mopts->krakOut) {
              writeAlignmentsToKrakenDump(read, formatter,
                                          validHits, bstream);
              alignmentStreamCount += validHits.size();
            } else if (mopts->salmonOut) {
              writeAlignmentsToKrakenDump(read, formatter,
                                            validHits, bstream, false);
              alignmentStreamCount += validHits.size();
            } else if (mopts->radOut) {
              writeAlignmentsToRADSingle(read, formatter,
                                            jointAlignments, bstream, false);
              alignmentStreamCount += validHits.size();
            } else if (jointHits.size() > 0 and !mopts->noOutput) {
                // write sam output for mapped reads
                writeAlignmentsToStreamSingle(read, formatter, jointAlignments, sstream, !mopts->noOrphan);
                alignmentStreamCount += jointAlignments.size();
            } else if (jointHits.size() == 0 and !mopts->noOutput) {
                // write sam output for un-mapped reads
              writeUnalignedSingleToStream(read, formatter, sstream);
              alignmentStreamCount += 1;
            }

            // write them on cmd
            if (hctr.numReads > hctr.lastPrint + 100000) {
                hctr.lastPrint.store(hctr.numReads.load());
                if (!mopts->quiet and iomutex->try_lock()) {
                    if (hctr.numReads > 0) {
                        std::cerr << "\r\r";
                    }
                    std::cerr << "saw " << hctr.numReads << " reads : "
                              << "pe / read = "
                              << hctr.peHits / static_cast<float>(hctr.numReads)
                              << " : se / read = "
                              << hctr.seHits / static_cast<float>(hctr.numReads) << ' ';
                    iomutex->unlock();
                }
            }

            // try dumping the output
            bool last_read = (read_it + 1 == rg.end());
            if (!mopts->noOutput and (alignmentStreamCount > alignmentStreamLimit or last_read)) {
                // Get rid of last newline
                if (mopts->krakOut || mopts->salmonOut) {
                    if (mopts->salmonOut && bstream.getBytes() > 0) {
                        BinWriter sbw(64);
                        sbw << bstream.getBytes();
                        outQueue->info("{}{}", sbw, bstream);
                    } else if (mopts->krakOut) {
                        outQueue->info("{}", bstream);
                    }
                    bstream.clear();
                } else if (mopts->radOut) {
                    if (mopts->salmonOut && bstream.getBytes() > 0) {
                        BinWriter sbw(64);
                        sbw << bstream.getBytes() << chunkReads;
                        outQueue->info("{}{}", sbw, bstream);
                        chunkReads = 0;
                    } else if (mopts->krakOut) {
                        outQueue->info("{}", bstream);
                    }
                    bstream.clear();
                } else {
                    std::string outStr(sstream.str());
                    // Get rid of last newline
                    if (!outStr.empty()) {
                        outStr.pop_back();
                        outQueue->info(std::move(outStr));
                    }
                    sstream.clear();
                }
                alignmentStreamCount = 0;
            }
        } // for all reads in this job
    } // processed all reads
}

//===========
// PAIRED END
//============
template<typename PufferfishIndexT>
bool spawnProcessReadsThreads(
        uint32_t nthread,
        paired_parser *parser,
        PufferfishIndexT &pfi,
        MutexT &iomutex,
        std::shared_ptr<spdlog::logger> outQueue,
        HitCounters &hctr,
        phmap::flat_hash_set<std::string>& gene_names,
        phmap::flat_hash_set<std::string>& rrna_names,
        pufferfish::AlignmentOpts *mopts) {

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nthread; ++i) {

        threads.emplace_back(processReadsPair<PufferfishIndexT>,
                             parser,
                             std::ref(pfi),
                             &iomutex,
                             outQueue,
                             std::ref(hctr),
                             std::ref(gene_names),
                             std::ref(rrna_names),
                             mopts);
    }
    for (auto &t : threads) { t.join(); }

    return true;
}

//===========
// SINGLE END
//============
template<typename PufferfishIndexT>
bool spawnProcessReadsThreads(
        uint32_t nthread,
        single_parser *parser,
        PufferfishIndexT &pfi,
        MutexT &iomutex,
        std::shared_ptr<spdlog::logger> outQueue,
        HitCounters &hctr,
        phmap::flat_hash_set<std::string>& gene_names,
        pufferfish::AlignmentOpts *mopts) {

    std::vector<std::thread> threads;

    for (size_t i = 0; i < nthread; ++i) {

        threads.emplace_back(processReadsSingle<PufferfishIndexT>,
                             parser,
                             std::ref(pfi),
                             &iomutex,
                             outQueue,
                             std::ref(hctr),
                             std::ref(gene_names),
                             mopts);
    }
    for (auto &t : threads) { t.join(); }

    return true;
}

void printAlignmentSummary(HitCounters &hctrs, std::shared_ptr<spdlog::logger> consoleLog) {
    consoleLog->info("Done mapping reads.");
    consoleLog->info("\n\n");
    consoleLog->info("=====");
    consoleLog->info("Observed {} reads", hctrs.numReads);
    consoleLog->info("Number of reads totally discarded for being smaller than k: {} read", hctrs.tooShortReads);
    consoleLog->info("Rate of Fragments with at least one found k-mer: {:03.2f}%",
                     (100.0 * static_cast<float>(hctrs.numMappedAtLeastAKmer)) / hctrs.numReads);
    consoleLog->info("Discordant Rate: {:03.2f}%",
                     (100.0 * static_cast<float>(hctrs.numOfOrphans)) / hctrs.numReads);
    consoleLog->info("Total reads Mapped: {}", (hctrs.numMapped));
    consoleLog->info("Mapping rate : {:03.2f}%", (100.0 * static_cast<float>(hctrs.numMapped)) / hctrs.numReads);
    consoleLog->info("Average # hits per read : {}", hctrs.totAlignment / static_cast<float>(hctrs.numReads));
    consoleLog->info("Total # of alignments : {}", hctrs.totAlignment);
    consoleLog->info("Total # of orphans : {}", hctrs.numOfOrphans);
    consoleLog->info("Total # of pe hits : {}", hctrs.peHits);
    consoleLog->info("Total # of total Hits : {}", hctrs.totHits);
    //consoleLog->info("Total # of valid hits : {}", hctrs.validHits);
    consoleLog->info("Max multimapping group : {}", hctrs.maxMultimapping);
    consoleLog->info("Total number of alignment attempts : {}", hctrs.totalAlignmentAttempts);
    consoleLog->info("Number of skipped alignments because of cache hits : {}", hctrs.skippedAlignments_byCache);
    consoleLog->info("Number of skipped alignments because of perfect chains : {}", hctrs.skippedAlignments_byCov);
    consoleLog->info("Number of alignments calculations skipped by non-alignable: {}", hctrs.skippedAlignments_notAlignable);

    consoleLog->info("Number of cigar strings which are fixed: {}", hctrs.cigar_fixed_count);
    consoleLog->info("=====");
}

template<typename PufferfishIndexT>
bool alignReads(
        PufferfishIndexT &pfi,
        std::shared_ptr<spdlog::logger> consoleLog,
        pufferfish::AlignmentOpts *mopts) {

    std::streambuf *outBuf;
    std::ofstream outFile;
    std::unique_ptr<std::ostream> outStream{nullptr};
    //bool haveOutputFile{false} ;
    std::shared_ptr<spdlog::logger> outLog{nullptr};

    phmap::flat_hash_set<std::string> gene_names;
    if (mopts->filterGenomics or mopts->filterMicrobiomBestScore) {
        std::ifstream gene_names_file;
        std::string gene_name;
        gene_names_file.open(mopts->genesNamesFile);
        if (!gene_names_file) {
            std::cerr << "Genomic file does not exist\n";
            exit(1);
        }
        std::string sepStr = " \t";
        while (gene_names_file >> gene_name) {
            auto processedName =
                    gene_name.substr(0, gene_name.find_first_of(sepStr));
            gene_names.insert(processedName);
        }
        gene_names_file.close();
    }

    phmap::flat_hash_set<std::string> rrna_names;
    if (mopts->filterMicrobiom) {
        std::ifstream rrna_names_file;
        std::string rrna_name;
        rrna_names_file.open(mopts->rrnaFile);
        if (!rrna_names_file) {
            std::cerr << "Genomic file does not exist\n";
            exit(1);
        }
        std::string sepStr = " \t";
        while (rrna_names_file >> rrna_name) {
            auto processedName =
                    rrna_name.substr(0, rrna_name.find_first_of(sepStr));
            rrna_names.insert(processedName);
        }
        rrna_names_file.close();
    }


    if (!mopts->noOutput) {
        if (mopts->outname == "-") {
            outBuf = std::cout.rdbuf();
        } else {
            outFile.open(mopts->outname);
            outBuf = outFile.rdbuf();
            //haveOutputFile = true ;
        }

        // out stream to the buffer
        // it can be std::cerr or a file
        if (mopts->compressedOutput) {
            outStream.reset(new zstr::ostream(outBuf));
        } else {
            outStream.reset(new std::ostream(outBuf));
        }
        // the async queue size must be a power of 2
        size_t queueSize{2*mopts->numThreads};
        spdlog::set_async_mode(queueSize);

        if (mopts->krakOut || mopts->salmonOut || mopts->radOut) {
            auto outputSink = std::make_shared<ostream_bin_sink_mt>(*outStream);
            outLog = std::make_shared<spdlog::logger>("puffer::outLog", outputSink);
            outLog->set_pattern("");
        } else {
            auto outputSink = std::make_shared<spdlog::sinks::ostream_sink_mt>(*outStream);
            outLog = std::make_shared<spdlog::logger>("puffer::outLog", outputSink);
            outLog->set_pattern("%v");
        }
        // write the SAM Header
        // If nothing gets printed by this time we are in troubleR
        if (mopts->radOut) {
            writeRADHeader(pfi, outLog, mopts);
        } else if (mopts->krakOut || mopts->salmonOut) {
            writeKrakOutHeader(pfi, outLog, mopts);
        } else { //TODO do we need to remove the txp from the list? The ids are then invalid
            writeSAMHeader(pfi, outLog,
                    mopts->filterGenomics or mopts->filterMicrobiom or mopts->filterMicrobiomBestScore,
                    gene_names,
                    rrna_names);
        }
    }

    uint32_t nthread = mopts->numThreads;
    std::unique_ptr<paired_parser> pairParserPtr{nullptr};
    std::unique_ptr<single_parser> singleParserPtr{nullptr};

    size_t chunkSize{10000};
    MutexT iomutex;

    if (!mopts->singleEnd) {
        ScopedTimer timer(!mopts->quiet);
        HitCounters hctrs;
        consoleLog->info("mapping reads ... \n\n\n");
        std::vector<std::string> read1Vec = pufferfish::util::tokenize(mopts->read1, ',');
        std::vector<std::string> read2Vec = pufferfish::util::tokenize(mopts->read2, ',');

        if (read1Vec.size() != read2Vec.size()) {
            consoleLog->error("The number of provided files for"
                              "-1 and -2 are not same!");
            std::exit(1);
        }

        uint32_t nprod = (read1Vec.size() > 1) ? 2 : 1;
        pairParserPtr.reset(new paired_parser(read1Vec, read2Vec, nthread, nprod, chunkSize));
        pairParserPtr->start();
        spawnProcessReadsThreads(nthread, pairParserPtr.get(), pfi, iomutex,
                                 outLog, hctrs, gene_names, rrna_names, mopts);
        pairParserPtr->stop();
        consoleLog->info("flushing output queue.");
        printAlignmentSummary(hctrs, consoleLog);
        if (outLog) { outLog->flush(); }
    } else {
        ScopedTimer timer(!mopts->quiet);
        HitCounters hctrs;
        consoleLog->info("mapping reads ... \n\n\n");
        std::vector<std::string> readVec = pufferfish::util::tokenize(mopts->unmatedReads, ',');

        uint32_t nprod = (readVec.size() > 1) ? 2 : 1;
        singleParserPtr.reset(new single_parser(readVec, nthread, nprod, chunkSize));
        singleParserPtr->start();

        spawnProcessReadsThreads(nthread, singleParserPtr.get(), pfi, iomutex,
                                 outLog, hctrs, gene_names, mopts);

        singleParserPtr->stop();
        consoleLog->info("flushing output queue.");
        printAlignmentSummary(hctrs, consoleLog);
        if (outLog) { outLog->flush(); }
    }
    return true;
}

template<typename PufferfishIndexT>
bool alignReadsWrapper(
        PufferfishIndexT &pfi,
        std::shared_ptr<spdlog::logger> consoleLog,
        pufferfish::AlignmentOpts *mopts) {
    bool res = true;
    if (mopts->listOfReads) {
        uint64_t readCntr = 1;
        if (mopts->singleEnd) {
            std::string unmatedReadsFile = mopts->unmatedReads;
            std::ifstream unmatedReadsF(unmatedReadsFile);
            std::string outname = mopts->outname;
            unmatedReadsF >> mopts->unmatedReads;
            while (unmatedReadsF.good() and mopts->unmatedReads != "") {
                consoleLog->info("Read {}: {}", readCntr, mopts->unmatedReads);
                readCntr++;
                uint64_t start = mopts->unmatedReads.find_last_of('/');
                if (start == std::string::npos) {
                    start = 0;
                } else {
                    start += 1;
                }
                uint64_t end = mopts->unmatedReads.find_last_of('.');
                mopts->outname = outname + mopts->unmatedReads.substr(start, end - start);
                if (mopts->salmonOut) {
                    mopts->outname += ".pam";
                } else {
                    mopts->outname += ".sam";
                }
                res &= alignReads(pfi, consoleLog, mopts);
                unmatedReadsF >> mopts->unmatedReads;
            }
        } else {
            std::string readFile1 = mopts->read1;
            std::string readFile2 = mopts->read2;
            std::ifstream readF1(readFile1);
            std::ifstream readF2(readFile2);
            std::string outname = mopts->outname;
            readF1 >> mopts->read1;
            readF2 >> mopts->read2;
            while (readF1.good() and readF2.good() and mopts->read1 != "" and mopts->read2 != "") {
                consoleLog->info("Read Pair {}: {}, {}", readCntr, mopts->read1, mopts->read2);
                readCntr++;
                uint64_t start = mopts->read1.find_last_of('/');
                if (start == std::string::npos) {
                    start = 0; // if / not found, start from index 0 of read name
                } else {
                    start+=1; // if / found, start the output file name from the index after /
                }
                uint64_t end = mopts->read1.find_last_of('.');
                if (end-2 == start) {
                    end-=1; // for the very rare names such as r1.fastq and r2.fastq
                } else {
                    end-=2; // otherwise, don't worry about the one last character missing for cases such as read1.fastq
                }
                mopts->outname = outname + mopts->read1.substr(start, end - start);
                if (mopts->salmonOut) {
                    mopts->outname += ".pam";
                } else {
                    mopts->outname += ".sam";
                }
                res &= alignReads(pfi, consoleLog, mopts);
                readF1 >> mopts->read1;
                readF2 >> mopts->read2;
            }
        }
    } else res &= alignReads(pfi, consoleLog, mopts);
    return res;
}

int pufferfishAligner(pufferfish::AlignmentOpts &alnargs) {

    auto consoleLog = spdlog::stderr_color_mt("console");
    bool success{false};
    auto indexDir = alnargs.indexDir;

    std::string indexType;
    {
        std::ifstream infoStream(indexDir + "/info.json");
        cereal::JSONInputArchive infoArchive(infoStream);
        infoArchive(cereal::make_nvp("sampling_type", indexType));
        std::cerr << "Index type = " << indexType << "\n";
        infoStream.close();
    }

    if (indexType == "dense") {
        PufferfishIndex pfi(indexDir);
        success = alignReadsWrapper(pfi, consoleLog, &alnargs);
    } else if (indexType == "sparse") {
        PufferfishSparseIndex pfi(indexDir);
        success = alignReadsWrapper(pfi, consoleLog, &alnargs);
    } else if (indexType == "lossy") {
        PufferfishLossyIndex pfi(indexDir);
        success = alignReadsWrapper(pfi, consoleLog, &alnargs);
    }

    if (!success) {
        consoleLog->warn("Problem mapping.");
    }
    return 0;
}
