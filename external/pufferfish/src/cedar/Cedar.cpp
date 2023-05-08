#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib> // std::abs
#include <cmath> // std::abs
#include <memory>
#include <unordered_set>
#include <list>

#include "SetCover.h"
#include "cedar/EquivalenceClassBuilder.hpp"
#include "cedar/Cedar.hpp"

#include "clipp.h"
#include "CLI/Timer.hpp"
#include "PufferFS.hpp"
#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/string.hpp"
//#include "LinearMultiArray.h"

#ifndef LIBASYNC_STATIC
    #define LIBASYNC_STATIC true
#endif

#include "async++.h"

struct set_data {
    unsigned int **sets;
    unsigned short **weights;
    unsigned int *set_sizes;
    unsigned int *element_size_lookup;
    unsigned int set_count;
    unsigned int uniqu_element_count;
    unsigned int all_element_count;
    unsigned int max_weight;
};

// trim from both ends (in place)
static inline std::string &trim(std::string &s) {
    std::string chars = "\t\n\v\f\r ";
    s.erase(0, s.find_first_not_of(chars));
    s.erase(s.find_last_not_of(chars) + 1);
    return s;
}

struct CedarOpts {
    std::string taxonomyTree_filename;
    std::string refId2TaxId_filename;
    std::string mapperOutput_filename;
    std::string output_filename;
    std::string level = "species";
    size_t maxIter = 1000;
    double eps = 0.001;
    double filterThreshold = 0;
    double minCnt = 0;
    bool flatAbund{false};
    bool requireConcordance{false};
    bool isPuffOut{false};
    bool isSAM{false};
    bool onlyUniq{false};
    bool onlyPerfect{false};
    uint32_t segmentSize{200};
    uint32_t rangeFactorizationBins{4};
    uint32_t numThreads{4};
};

template<class ReaderType, class FileReaderType>
Cedar<ReaderType, FileReaderType>::Cedar(std::string &taxonomyTree_filename,
                         std::string &refId2TaxId_filename,
                         std::string pruneLevelIn,
                         double filteringThresholdIn,
                         bool flatAbundIn,
                         std::shared_ptr<spdlog::logger> loggerIn) {
    flatAbund = flatAbundIn;
    logger = loggerIn;
    logger->info("Cedar: Construct ..");

    // map rank string values to enum values
    filteringThreshold = filteringThresholdIn;
    if (!flatAbund) {
        pruningLevel = TaxaNode::str2rank(pruneLevelIn);
        std::ifstream tfile;
        uint32_t id, pid;
        std::string rank, name;
        // load the reference id (name) to its corresponding taxonomy id
        tfile.open(refId2TaxId_filename);
        while (!tfile.eof()) {
            tfile >> name >> id;
            refId2taxId[name] = id;
        }
        tfile.close();

        // load the taxonomy child-parent tree and the rank of each node
        tfile.open(taxonomyTree_filename);
        std::string tmp, line;
        uint64_t lcntr{0};
        //std::cerr << "\n\nreading taxonomy tree\n\n";
        while (!tfile.eof()) {
            try {
                lcntr++;
                std::getline(tfile, line);
                uint64_t first = line.find_first_of('|');
                std::string fstr = line.substr(0, first);
                id = (uint32_t) std::stoul(trim(fstr));
                uint64_t second = line.find_first_of('|', first + 1);
                fstr = line.substr(first + 1, second - first - 1);
                pid = (uint32_t) std::stoul(trim(fstr));
                uint64_t third = line.find_first_of('|', second + 1);
                fstr = line.substr(second + 1, third - second - 1);
                rank = trim(fstr);
                if (rank.find("HS") != std::string::npos)
                    std::cerr << "here:" << id << " " << pid << " " << tmp << "\n";
                taxaNodeMap[id] = TaxaNode(id, pid, TaxaNode::str2rank(rank));
                if (taxaNodeMap[id].isRoot()) {
                    rootId = id;
                    logger->info("Root Id : {}", id);
                }
            } catch (const std::invalid_argument &ia) {
                std::cerr << "Invalid argument: " << ia.what() << '\n';
                continue;
            }

        }

        tfile.close();
    }
}

template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::calculateCoverage() {
    for (uint32_t i = 0; i < strain_coverage_bins.size(); i++) {
        auto bins = strain_coverage_bins[i];
        double covered = 0.0;
        uint32_t expression_depth = 0;
        for (uint32_t j = 0; j < bins.size(); j++) {
            if (bins[j] > 0)
                covered++;
            expression_depth += bins[j];
        }
        strain_coverage[i] = covered / bins.size();
    }
}

template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::processAlignmentBatch(uint32_t threadID,
                                              std::vector<ReadInfo> &alignmentGrp,
                                              Stats &stats,
                                              EquivalenceClassBuilder& eqb,
                                              std::mutex &iomutex,
                                              bool requireConcordance,
                                              bool onlyUniq,
                                              bool onlyPerfect,
                                              uint32_t segmentSize,
                                              uint32_t rangeFactorizationBins,
                                              bool getReadName) {
    TaxaNode *prevTaxa{nullptr};
    std::set<uint64_t> activeTaxa;
    ReaderType mappings(fileReader.getReader(), logger);
    while (mappings.nextAlignmentGroup(alignmentGrp, iomutex, threadID, getReadName)) {
        for (auto &readInfo: alignmentGrp) {
            stats.totalReadCnt++;
            if (stats.totalReadCnt % 1000000 == 0) {
//            {
//                std::lock_guard<std::mutex> l(iomutex);
                std::cerr << "\nThread " << threadID << " processed " << stats.totalReadCnt << " reads";
//            }
            }
            activeTaxa.clear();
            double readMappingsScoreSum = 0;
            std::vector<std::pair<uint64_t, double>> readPerStrainProbInst;
            readPerStrainProbInst.reserve(readInfo.cnt);
            bool isConflicting = true;
            uint64_t maxScore = readInfo.len;
            //std::cerr << "Thread " << threadID << " " << stats.totalReadCnt << " " <<readInfo.cnt << "\n";
            if (readInfo.cnt != 0) {
                if (readInfo.cnt > 1) {
                    stats.totalMultiMappedReads++;
                }

                std::set<uint64_t> seen;
                prevTaxa = nullptr;
                for (auto &mapping : readInfo.mappings) {
                    auto &refNam = fileReader.refName(mapping.getId());
                    // first condition: Ignore those references that we don't have a taxaId
                    // second condition: Ignore repeated exactly identical mappings (FIXME thing)
                    // note: Fatemeh!! don't change the or and and again!!
                    // if we are on flatAbund, we want to count for multiple occurrences of a reference
                    if ((flatAbund or
                         (refId2taxId.find(refNam) != refId2taxId.end()))
                        and
                        activeTaxa.find(mapping.getId()) == activeTaxa.end()) {
                        if (prevTaxa != nullptr and
                            activeTaxa.find(mapping.getId()) == activeTaxa.end() and
                            !prevTaxa->compareIntervals(mapping)) {
                            isConflicting = false;
                        }
                        activeTaxa.insert(mapping.getId());
                        if (requireConcordance && fileReader.isMappingPaired() &&
                            (!mapping.isConcordant() || mapping.isFw(ReadEnd::LEFT) == mapping.isFw(ReadEnd::RIGHT))) {
                            stats.discordantMappings++;
                            continue;
                        }

                        /*auto tid = flatAbund ? mapping.getId() : refId2taxId[refNam];
                        seqToTaxMap[mapping.getId()] = static_cast<uint32_t>(tid);*/

                       /* if (cov.find(refId2taxId[refNam]) == cov.end()) {
                            util::update(cov[refId2taxId[refNam]], 0);
                        }*/
                        auto tid = flatAbund ? mapping.getId() : refId2taxId[refNam];
                        //std::cerr << "c0:" << cov[tid] << " " << mapping.getScore() << " ";
                        if (cov.contains(tid) == 0) {
                            std::cerr << tid <<" doesn't exist in cov\n";
                            std::exit(1);
                        }
                        util::incLoop(cov[tid], mapping.getScore());
                        //std::cerr << "c1:" << cov[tid] << "\n";
                        readPerStrainProbInst.emplace_back(mapping.getId(),
                                                           static_cast<double>( mapping.getScore()) /
                                                           static_cast<double>(fileReader.refLength(mapping.getId())));
                        readMappingsScoreSum += readPerStrainProbInst.back().second;

                        uint32_t bin_number = mapping.getPos(ReadEnd::LEFT) / segmentSize > 1 ?
                                              mapping.getPos(ReadEnd::LEFT) / segmentSize - 1 : 0;
                        auto &cbins = strain_coverage_bins[mapping.getId()];
                        if (bin_number < cbins.size()) {
                            cbins[bin_number]++;
                        } else {
                            std::cerr << "SHOULDN'T HAPPEN! out of bound: " << mapping.getId() << " " << refNam << " " << fileReader.refLength(mapping.getId()) << " " << bin_number << " " << cbins.size() << " ";
                            std::cerr << "lpos: " << mapping.getPos(ReadEnd::LEFT) << "\n";
                        }
                        prevTaxa = &mapping;
                    }
                }

                if (activeTaxa.empty()) {
                    stats.seqNotFound++;
                } else if ((!onlyUniq and !onlyPerfect)
                           or (onlyUniq and activeTaxa.size() == 1)
                           or (onlyPerfect and activeTaxa.size() == 1
                               and
                               readInfo.mappings[0].getScore() >= maxScore)) {
                    if (!isConflicting) { stats.conflicting++; }
                    // it->first : strain id
                    // it->second : prob of current read comming from this strain id (mapping_score/ref_len)
                    double probsum{0.0};
                    for (auto it = readPerStrainProbInst.begin(); it != readPerStrainProbInst.end(); it++) {
                        it->second = it->second / readMappingsScoreSum; // normalize the probabilities for each read
                        if (strain.contains(it->first) == 0) {
                            std::cerr << it->first << " doesn't exist in strain\n";
                            std::exit(1);
                        }
                        util::incLoop(strain[it->first], 1.0 / static_cast<double>(readPerStrainProbInst.size()));
                        probsum += 1.0 / static_cast<double>(readPerStrainProbInst.size());
                    }
                    stats.globalprobsum += probsum;
                    if (abs(probsum - 1.0) >= 1e-10) {
                        std::cerr << "hfs!!";
                        std::exit(1);
                    }
                    // SAVE MEMORY, don't push this
                    //readPerStrainProb.push_back(readPerStrainProbInst);

                    // construct the range factorized eq class here
                    std::sort(readPerStrainProbInst.begin(), readPerStrainProbInst.end(),
                              [](std::pair<uint64_t, double> &a, std::pair<uint64_t, double> &b) {
                                  return a.first < b.first;
                              });
                    std::vector<uint32_t> genomeIDs;
                    genomeIDs.reserve(2 * readPerStrainProbInst.size());
                    std::vector<double> probs;
                    probs.reserve(readPerStrainProbInst.size());
                    for (auto &it : readPerStrainProbInst) {
                        genomeIDs.push_back(static_cast<uint32_t>(it.first));
                        probs.push_back(it.second);
                    }
                    if (rangeFactorizationBins > 0) {
                        uint64_t genomeSize = genomeIDs.size();
                        uint64_t rangeCount =
                                static_cast<uint64_t>(std::sqrt(genomeSize)) + rangeFactorizationBins;
                        for (uint64_t i = 0; i < genomeSize; i++) {
                            int rangeNumber = static_cast<int>(probs[i] * rangeCount);
                            genomeIDs.push_back(static_cast<uint32_t>(rangeNumber));
                        }
                    }
                    stats.readCnt++;
                    TargetGroup tg(genomeIDs);
                    eqb.addGroup(std::move(tg), probs); //add or update eq read cnt by 1
                } else {
                    stats.totalReadsNotPassingCond++;
                }
            } else {
                stats.totalUnmappedReads++;
            }
        }
    }
}


template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::loadMappingInfo(std::string mapperOutput_filename,
                                        bool requireConcordance,
                                        bool onlyUniq,
                                        bool onlyPerfect,
                                        uint32_t segmentSize,
                                        uint32_t rangeFactorizationBins,
                                        uint32_t nThreads) {

    Stats stats;
    logger->info("Cedar: Load Mapping File ..");
    logger->info("Mapping Output File: {}", mapperOutput_filename);
    fileReader.load(mapperOutput_filename, logger);
    logger->info("is dataset paired end? {}", fileReader.isMappingPaired());

    // Construct coverage bins
    for (uint64_t i = 0; i < fileReader.numRefs(); ++i) {
        auto refLength = fileReader.refLength(i);
        uint64_t binCnt = refLength / segmentSize;
        if (binCnt == 0) binCnt = 1;
        std::vector<uint32_t> bins(binCnt, 0);
        strain_coverage_bins[i] = bins;
        strain[i] = 0;
        auto tid = flatAbund ? i : refId2taxId[fileReader.refName(i)];
        seqToTaxMap[i] = static_cast<uint32_t>(tid);
        cov[tid] = 0;
    }

    constexpr const bool getReadName = true;
    std::cerr << "\n";


    std::vector<std::thread> threads;
    std::vector<Stats> statsPerThread(nThreads);
    std::vector<EquivalenceClassBuilder> eqbPerThread(nThreads);
    std::vector<std::vector<ReadInfo>> alignmentGroups(nThreads);
    for (uint32_t i = 0; i < nThreads; ++i) alignmentGroups[i].resize(ALIGNMENTS_PER_BATCH);
    std::mutex iomutex;
    for (uint32_t i = 0; i < nThreads; ++i) {
            auto threadFun = [&, i]() -> void {
            processAlignmentBatch(i,
                                  alignmentGroups[i],
                                  statsPerThread[i],
                                  eqbPerThread[i],
                                  iomutex,
                                  requireConcordance,
                                  onlyUniq,
                                  onlyPerfect,
                                  segmentSize,
                                  rangeFactorizationBins,
                                  getReadName);
        };
        threads.emplace_back(threadFun);
    }
    for (auto &t : threads) {
        t.join();
    }
    for (auto& s: statsPerThread) {stats.update(s);}
    readCnt = stats.readCnt;
   /* for (auto &s: strains) {
        for (auto& kv: s) {
            if (strain.find(kv.first) == strain.end())
                strain[kv.first] = 0;
            strain[kv.first] += kv.second;
        }
    }
    for (auto &s: covs) {
        for (auto& kv: s) {
            if (cov.find(kv.first) == cov.end())
                cov[kv.first] = 0;
            cov[kv.first] += kv.second;
        }
    }*/
    for (auto &eq: eqbPerThread) {
        eqb.mergeUnfinishedEQB(eq);
    }
    calculateCoverage();
    std::cerr << "\r";
    //logger->info("Total # of unique reads: {}", readset.size());
    //notMappedReadsFile.close();
    logger->info("# of mapped (and accepted) reads: {}", stats.readCnt);
    logger->info("global probsum: {}", stats.globalprobsum);
    if (onlyUniq or onlyPerfect)
        logger->info("# of mapped reads that were not uniq/perfect: {}", stats.totalReadsNotPassingCond);
    logger->info("# of multi-mapped reads: {}", stats.totalMultiMappedReads);
    logger->info("# of conflicting reads: {}", stats.conflicting);
    logger->info("# of unmapped reads: {}", stats.totalUnmappedReads);
    if (requireConcordance)
        logger->info("Discarded {} discordant mappings.", stats.discordantMappings);
}

template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::loadMappingInfo(std::string mapperOutput_filename,
                                        bool requireConcordance,
                                        bool onlyUniq,
                                        bool onlyPerfect,
                                        uint32_t segmentSize,
                                        uint32_t rangeFactorizationBins) {
    uint32_t rangeFactorization{rangeFactorizationBins};
    uint64_t totalReadCnt{0}, seqNotFound{0},
            totalMultiMappedReads{0}, totalUnmappedReads{0}, totalReadsNotPassingCond{0}, tid;
    logger->info("Cedar: Load Mapping File ..");
    logger->info("Mapping Output File: {}", mapperOutput_filename);
    fileReader.load(mapperOutput_filename, logger);
    logger->info("is dataset paired end? {}", fileReader.isMappingPaired());
    ReadInfo readInfo;
    TaxaNode *prevTaxa{nullptr};
    size_t conflicting{0};
    size_t discordantMappings{0};

    // Construct coverage bins
    for (uint64_t i = 0; i < fileReader.numRefs(); ++i) {
        auto refLength = fileReader.refLength(i);
        uint64_t binCnt = refLength / segmentSize;
        if (binCnt == 0) binCnt = 1;
        std::vector<uint32_t> bins(binCnt, 0);
        strain_coverage_bins[i] = bins;
    }

    constexpr const bool getReadName = true;
    double globalprobsum{0.0};
    std::cerr << "\n";
    std::mutex iomutex;
    while (fileReader.nextRead(readInfo, iomutex, getReadName)) {
        totalReadCnt++;
        if (totalReadCnt % 1000000 == 0) {
            std::cerr << "\rProcessed " << totalReadCnt << " reads";
        }
        activeTaxa.clear();
        double readMappingsScoreSum = 0;
        std::vector<std::pair<uint64_t, double>> readPerStrainProbInst;
        readPerStrainProbInst.reserve(readInfo.cnt);
        bool isConflicting = true;
        uint64_t maxScore = readInfo.len;
        if (readInfo.cnt != 0) {
            if (readInfo.cnt > 1) {
                totalMultiMappedReads++;
            }

            std::set<uint64_t> seen;
            prevTaxa = nullptr;
            std::vector<std::tuple<uint32_t, bool> > mapping_scores;
            for (auto &mapping : readInfo.mappings) {
                auto &refNam = fileReader.refName(mapping.getId());
                // first condition: Ignore those references that we don't have a taxaId
                // second condition: Ignore repeated exactly identical mappings (FIXME thing)
                // note: Fatemeh!! don't change the or and and again!!
                // if we are on flatAbund, we want to count for multiple occurrences of a reference
                if ((flatAbund or
                     (refId2taxId.find(refNam) != refId2taxId.end()))
                    and
                    activeTaxa.find(mapping.getId()) == activeTaxa.end()) {
                    if (prevTaxa != nullptr and
                        activeTaxa.find(mapping.getId()) == activeTaxa.end() and
                        !prevTaxa->compareIntervals(mapping)) {
                        isConflicting = false;
                    }
                    activeTaxa.insert(mapping.getId());
                    if (requireConcordance && fileReader.isMappingPaired() &&
                        (!mapping.isConcordant() || mapping.isFw(ReadEnd::LEFT) == mapping.isFw(ReadEnd::RIGHT))) {
                        discordantMappings++;
                        continue;
                    }

                    tid = flatAbund ? mapping.getId() : refId2taxId[refNam];
                    seqToTaxMap[mapping.getId()] = static_cast<uint32_t >(tid);

                    if (cov.find(refId2taxId[refNam]) == cov.end()) {
                        cov[refId2taxId[refNam]] = 0;
                    }
                    mapping_scores.push_back(std::tuple<uint32_t, bool>(mapping.getScore(), mapping.isPaired()));
                    util::incLoop(cov[refId2taxId[refNam]], mapping.getScore());
                    readPerStrainProbInst.emplace_back(mapping.getId(),
                                                       static_cast<double>( mapping.getScore()) /
                                                       static_cast<double>(fileReader.refLength(mapping.getId())));
                    readMappingsScoreSum += readPerStrainProbInst.back().second;

                    uint32_t bin_number = mapping.getPos(ReadEnd::LEFT) / segmentSize > 0 ?
                                          mapping.getPos(ReadEnd::LEFT) / segmentSize - 1 : 0;
                    auto &cbins = strain_coverage_bins[mapping.getId()];
                    if (bin_number < cbins.size()) {
                        cbins[bin_number]++;
                    } else {
                        std::cerr << "SHOULDN'T HAPPEN! out of bound: " << bin_number << " " << cbins.size() << " ";
                        std::cerr << "lpos: " << mapping.getPos(ReadEnd::LEFT) << "\n";
                    }
                    prevTaxa = &mapping;
                }
            }
            /*if (mapping_scores.size() > 20) {
                std::sort(mapping_scores.begin(),mapping_scores.end());
                std::cout<< readInfo.len << "\n";
                for (auto &mapping_score : mapping_scores) {
                    std::cout << std::get<0>(mapping_score) << " " << std::get<1>(mapping_score) << " ";
                }
                std::cout<<"\n";
            }*/

            if (activeTaxa.empty()) {
                seqNotFound++;
            } else if ((!onlyUniq and !onlyPerfect)
                       or (onlyUniq and activeTaxa.size() == 1)
                       or (onlyPerfect and activeTaxa.size() == 1
                           and
                           readInfo.mappings[0].getScore() >= maxScore)) {
                if (!isConflicting) { conflicting++; }
                // it->first : strain id
                // it->second : prob of current read comming from this strain id (mapping_score/ref_len)
                double probsum{0.0};
                for (auto it = readPerStrainProbInst.begin(); it != readPerStrainProbInst.end(); it++) {
                    it->second = it->second / readMappingsScoreSum; // normalize the probabilities for each read
                    util::incLoop(strain[it->first], 1.0 / static_cast<double>(readPerStrainProbInst.size()));
                    probsum += 1.0 / static_cast<double>(readPerStrainProbInst.size());
                }
                globalprobsum += probsum;
                if (abs(probsum - 1.0) >= 1e-10) {
                    std::cerr << "hfs!!";
                    std::exit(1);
                }
                // SAVE MEMORY, don't push this
                //readPerStrainProb.push_back(readPerStrainProbInst);

                // construct the range factorized eq class here 
                std::sort(readPerStrainProbInst.begin(), readPerStrainProbInst.end(),
                          [](std::pair<uint64_t, double> &a, std::pair<uint64_t, double> &b) {
                              return a.first < b.first;
                          });
                std::vector<uint32_t> genomeIDs;
                genomeIDs.reserve(2 * readPerStrainProbInst.size());
                std::vector<double> probs;
                probs.reserve(readPerStrainProbInst.size());
                for (auto &it : readPerStrainProbInst) {
                    genomeIDs.push_back(static_cast<uint32_t>(it.first));
                    probs.push_back(it.second);
                }
                if (rangeFactorization > 0) {
                    int genomeSize = genomeIDs.size();
                    int rangeCount = std::sqrt(genomeSize) + rangeFactorization;
                    for (int i = 0; i < genomeSize; i++) {
                        int rangeNumber = static_cast<int>(probs[i] * rangeCount);
                        genomeIDs.push_back(static_cast<uint32_t>(rangeNumber));
                    }
                }
                readCnt++;
                TargetGroup tg(genomeIDs);
                eqb.addGroup(std::move(tg), probs); //add or update eq read cnt by 1
            } else {
                totalReadsNotPassingCond++;
            }
        } else {
            totalUnmappedReads++;
        }
    }
    calculateCoverage();
    std::cerr << "\r";
    //logger->info("Total # of unique reads: {}", readset.size());
    //notMappedReadsFile.close();
    logger->info("# of mapped (and accepted) reads: {}", readCnt);
    logger->info("global probsum: {}", globalprobsum);
    if (onlyUniq or onlyPerfect)
        logger->info("# of mapped reads that were not uniq/perfect: {}", totalReadsNotPassingCond);
    logger->info("# of multi-mapped reads: {}", totalMultiMappedReads);
    logger->info("# of conflicting reads: {}", conflicting);
    logger->info("# of unmapped reads: {}", totalUnmappedReads);
    if (requireConcordance)
        logger->info("Discarded {} discordant mappings.", discordantMappings);
}

template<class ReaderType, class FileReaderType>
bool Cedar<ReaderType, FileReaderType>::applySetCover(std::vector<double> &strainCnt,
                                      std::vector<bool> &strainValid,
                                      std::vector<bool> &strainPotentiallyRemovable,
                                      double minCnt,
                                      bool canHelp,
                                      bool verbose){
    uint64_t previouslyValid{0};
    auto &eqvec = eqb.eqVec();
    // rule to find potentially removable strains
    for (size_t i = 0; i < strainCnt.size(); ++i) {
        strainPotentiallyRemovable[i] = strainValid[i] ? strainCnt[i] <= minCnt : false;
        previouslyValid += strainValid[i];
    }
    // find list of all references with a unique *valid* read mapped to them
    std::unordered_set<uint64_t> uniqueReadRefs;
    for (auto &eqc : eqvec) {
        auto &tg = eqc.first;
        auto &v = eqc.second;
        auto csize = v.weights.size();
        std::vector<uint64_t> eqValids;
        for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
            auto &tgt = tg.tgts[readMappingCntr];
            if (strainValid[tgt])
                eqValids.push_back(tgt);
        }
        if (eqValids.size() == 1) {
            uniqueReadRefs.insert(eqValids[0]);
        }
    }
    for (size_t i = 0; i < strainCnt.size(); ++i) {
        if (strainValid[i] and uniqueReadRefs.find(i) == uniqueReadRefs.end())
            strainPotentiallyRemovable[i] = true;
    }


    std::/*unordered_*/map<uint32_t, std::/*unordered_*/set<uint64_t>> ref2eqset;
    for (auto &eqc : eqvec) {
        auto &tg = eqc.first;
        auto &v = eqc.second;
        auto csize = v.weights.size();
        uint64_t totalValidRefsInCurEq{0}, potentialRemovableCntr{0};

        // count total number of valid references for each eq.
        // (other than those that have been invalidated in previous rounds)
        for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
            auto &tgt = tg.tgts[readMappingCntr];
            if (strainValid[tgt]) {
                totalValidRefsInCurEq++;
            }
            if (strainPotentiallyRemovable[tgt]) {
                potentialRemovableCntr++;
            }
        }
        if (potentialRemovableCntr >=
            totalValidRefsInCurEq) { // if all the refs in the eq are set as potentiallyRemovable
            for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
                auto &tgt = tg.tgts[readMappingCntr];
                if (strainValid[tgt] and strainPotentiallyRemovable[tgt]) {
                    ref2eqset[tgt].insert(tg.hash);
                }
            }
        }
    }

    std::unordered_map<uint64_t, uint32_t> eq2id;
// ref2eqset is a map from a reference to the set of ambiguous eqs it belongs to (careful, not all the set of eqs.)
// ambiguous eq is an eq that all its remaining refs has been chosen as potentially removable
    if (ref2eqset.size() > 0) {
// set cover input preparation
// convert the input to proper format for the library that runs setCover algo.
        uint32_t id = 1;
        for (auto &kv: ref2eqset) {
            for (auto v: kv.second)
                if (eq2id.find(v) == eq2id.end()) {
                    eq2id[v] = id;
                    id++;
                }
        }
        auto elementCnt = eq2id.size(); // n
        auto setCnt = ref2eqset.size(); // m
        set_data ret_struct;
        ret_struct.set_count = setCnt;
        ret_struct.uniqu_element_count = elementCnt;

        unsigned int *element_size = new unsigned int[elementCnt + 2];
        memset(&element_size[0], 0, sizeof(unsigned int) * (elementCnt + 2));
        ret_struct.
                element_size_lookup = element_size;
        unsigned int *set_size = new unsigned int[setCnt];

        ret_struct.set_sizes = set_size;
        ret_struct.sets = new unsigned int *[setCnt];
        ret_struct.weights = new unsigned short *[setCnt];
        ret_struct.max_weight = 0;
        ret_struct.all_element_count = 0;

        uint32_t i = 0;
        std::vector<uint64_t> setCoverId2ref(ref2eqset.size());
        std::vector<uint32_t> setWeight;
        for (auto &kv: ref2eqset) {
            setCoverId2ref[i] = kv.first;
            auto setSize = kv.second.size();
            ret_struct.set_sizes[i] = setSize;
            auto mw = static_cast<uint32_t>(this->strain_coverage[i] * 10000);
            auto perElement = static_cast<uint32_t>(mw/setSize);
            if (perElement < 1) {
                perElement = 1;
            }
            mw = perElement * setSize;
            setWeight.push_back(mw);
//            std::cerr << "pE: " << perElement << " " << mw << "\n";
            if (ret_struct.max_weight < mw/*setSize*/) {
                ret_struct.max_weight = mw;//setSize;
            }
            ret_struct.sets[i] = new unsigned int[setSize];
            uint32_t element_counter = 0;
            for (auto &eq: kv.second) {
                ret_struct.sets[i][element_counter++] = eq2id[eq];
                ret_struct.element_size_lookup[eq2id[eq]]++;
                ret_struct.all_element_count++;
            }
            ret_struct.weights[i] = new unsigned short[setSize];
            std::fill_n(ret_struct.weights[i], setSize, perElement);
            i++;
        }
// end of set cover input preparation
        if (verbose)
            std::cerr << "Set Cover Input:\n# of refs: " << ret_struct.set_count
                      << ", # of uniq eqs:" << ret_struct.uniqu_element_count
                      << ", max set size: " << ret_struct.max_weight << "\n";

// run setCover algo over the list of refs and eqs in ref2eqset
        set_cover setcover(ret_struct.set_count,
                           ret_struct.uniqu_element_count,
                           ret_struct.max_weight,
                           ret_struct.all_element_count,
                           ret_struct.element_size_lookup
        );
        if (verbose)
            std::cerr << "Done with executing setCover\n";
        for (uint64_t j = 0; j < ret_struct.set_count; j++) {
            setcover.add_set(j + 1, setWeight[j]/*ret_struct.set_sizes[j]*/,
                    (const unsigned int *) ret_struct.sets[j],
                             (const unsigned short *) ret_struct.weights[j],
                             ret_struct.set_sizes[j]);
            delete[] ret_struct.sets[j];
            delete[] ret_struct.weights[j];
//            free(ret_struct.sets[j]);
//            free(ret_struct.weights[j]);
        }
        if (verbose)
            std::cerr << "Done with executing setCover\n";
        std::list<set *> ret = setcover.execute_set_cover();
// put the list of minimum # of references that can cover all eqs (output of setCover algo.) in remainingRefs
        std::unordered_set<uint32_t> remainingRefs;
        for (auto iterator = ret.begin(); iterator != ret.end(); ++iterator) {
            remainingRefs.insert(setCoverId2ref[(*iterator)->set_id - 1]);
        }
        delete[] element_size;
        delete[] set_size;
        delete[] ret_struct.sets;
        delete[] ret_struct.weights;
// go over the list of references
        i = 0;
        for (uint64_t refCntr = 0; refCntr < strainValid.size(); refCntr++) {
            if (strainPotentiallyRemovable[refCntr] and
                remainingRefs.find(static_cast<unsigned int>(refCntr)) == remainingRefs.end()) {
                strainValid[refCntr] = false;
            }
            strainPotentiallyRemovable[refCntr] = false;
        }
        uint32_t totalValid{0};
        for (auto s: strainValid) {
            totalValid += s;
        }
        if (verbose)
            std::cerr << "Input refs: " << ref2eqset.size()
                      << "\tRescued refs: " << remainingRefs.size()
                      << "\n"
                      << "valid: " << totalValid
                      << ", invalid: " << strainValid.size() - totalValid << "\n\n";
        canHelp = previouslyValid != totalValid;
    }
    return canHelp;
}

template<class ReaderType, class FileReaderType>
bool Cedar<ReaderType, FileReaderType>::basicEM(size_t maxIter, double eps, double minCnt, uint32_t numThreads, bool verbose) {
    eqb.finish();
    auto &eqvec = eqb.eqVec();
    int64_t maxSeqID{-1};
    for (auto &kv : strain) {
        maxSeqID = (static_cast<int64_t>(kv.first) > maxSeqID) ? static_cast<int64_t>(kv.first) : maxSeqID;
    }

    std::vector<std::atomic < double>> newStrainCnt(maxSeqID + 1);
    for (auto i = 0; i < maxSeqID+1; i++) {
        util::update(newStrainCnt[i], 0.0);
    }
    std::vector<double> strainCnt(maxSeqID + 1);
    std::vector<bool> strainValid(maxSeqID + 1, true);
    std::vector<bool> strainPotentiallyRemovable(maxSeqID + 1, false);

    double validTaxIdCnt{0};
    for (auto &kv : strain) {
        strainCnt[kv.first] = kv.second;
        validTaxIdCnt += kv.second;
        if (fileReader.refName(kv.first) == "ENST00000357898.3|ENSG00000171621.13|OTTHUMG00000001279.4|OTTHUMT00000303969.1|SPSB1-202|SPSB1|1439|protein_coding|" or
            fileReader.refName(kv.first) == "ENST00000328089.10|ENSG00000171621.13|OTTHUMG00000001279.4|OTTHUMT00000003727.2|SPSB1-201|SPSB1|3120|protein_coding|") {
            std::stringstream ss;
            ss << kv.first << " " << fileReader.refName(kv.first) << " " << strainCnt[kv.first] << "\n";
            std::cerr << ss.str();
        }
    }
    logger->info("maxSeqID : {}", maxSeqID);
    logger->info("found : {} equivalence classes", eqvec.size());

    size_t totCount{0};
    for (auto &eqc : eqvec) {
        totCount += eqc.second.count;
    }
    logger->info("Total starting count {}", totCount);
    logger->info("Total mapped reads cnt {}", readCnt);
    logger->info("Total reads cnt mapped to valid taxids {}", static_cast<uint64_t >(validTaxIdCnt));

    size_t cntr = 0;
    bool converged = false;
    uint64_t thresholdingIterStep = 10;
    bool canHelp = true;
    //tbb::task_scheduler_init asyncScheduler(numThreads);
    async::threadpool_scheduler asyncScheduler(numThreads);
    uint64_t blockSize = 1000;
    std::vector<std::pair<decltype(eqvec.begin()), decltype(eqvec.begin())>>
        range(std::ceil(static_cast<double >(eqvec.size())/ static_cast<double>(blockSize)));
    uint64_t idx{0};
    for (uint64_t i = 0; i < eqvec.size(); i+=blockSize) {
        decltype(eqvec.begin()) rangeBegin(eqvec.begin()+i);
        auto rangeEnd = rangeBegin+blockSize > eqvec.end()? eqvec.end():rangeBegin+blockSize;
        range[idx] = std::make_pair(rangeBegin, rangeEnd);
        idx++;
    }
    while (cntr++ < maxIter && !converged) {
        if (cntr % thresholdingIterStep == 0 && canHelp) {
            canHelp = applySetCover(strainCnt, strainValid, strainPotentiallyRemovable, minCnt, canHelp, verbose);
        }
        // M step
        // Find the best (most likely) count assignment
        async::parallel_for(asyncScheduler,
                //tbb::blocked_range<size_t>(0, eqvec.size()),
//                            async::make_range(range, range+eqvec.size()),
                async::make_range(range.begin(), range.end()),
                [&eqvec, &strainValid, &strainCnt, &newStrainCnt, this]
                        (decltype(*range.begin())& rangePair/*const async::range<size_t> &range*/) -> void {
                    for (auto eqIt = rangePair.first; eqIt != rangePair.second; ++eqIt) {
                        auto &eqc = (*eqIt);//eqvec[eqID];
//                        std::cerr << eqID << "\n";
                        auto &tg = eqc.first;
                        auto &v = eqc.second;
                        auto csize = v.weights.size();
                        std::vector<double> tmpReadProb(csize, 0.0);
                        double denom{0.0};
                        for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
                            auto &tgt = tg.tgts[readMappingCntr];
                            if (strainValid[tgt]) {

                                tmpReadProb[readMappingCntr] =
                                        v.weights[readMappingCntr] * strainCnt[tgt] * this->strain_coverage[tgt];
                                // * (1.0/refLengths[tgt]);
                                denom += tmpReadProb[readMappingCntr];
                            }
                        }
                        for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
                            auto &tgt = tg.tgts[readMappingCntr];
                            if (strainValid[tgt]) {
//                                std::cerr << v.count * (tmpReadProb[readMappingCntr] / denom) << "\n";
                                util::incLoop(newStrainCnt[tgt], v.count * (tmpReadProb[readMappingCntr] / denom));
                            }
                            /*else {
                                newStrainCnt[tgt] = 0;
                            }*/
                        }
                    }
                });

        // E step
        // normalize strain probabilities using the denum : p(s) = (count(s)/total_read_cnt) 
        double readCntValidator = 0;
        converged = true;
        double maxDiff{0.0};
        // It's a small loop with critical variables so better not to be parallelized
        // It was a few seconds faster without parallelizing
        for (uint64_t i = 0; i < newStrainCnt.size(); ++i) {
            if (cntr < 20)
                if (fileReader.refName(i) == "ENST00000495701.1|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000099595.1|CEP104-211|CEP104|485|retained_intron|" or
                    fileReader.refName(i) == "ENST00000460038.5|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000009813.1|CEP104-206|CEP104|696|processed_transcript|" or
                    fileReader.refName(i) == "ENST00000494653.5|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000009749.1|CEP104-209|CEP104|2174|retained_intron|" or
                    fileReader.refName(i) == "ENST00000378230.7|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000009747.3|CEP104-202|CEP104|6424|protein_coding|") {
                std::stringstream ss;
                ss << i << " " << strainCnt[i] << "\n";
                std::cerr << ss.str();
            }
            readCntValidator += newStrainCnt[i];
            auto adiff = std::abs(newStrainCnt[i] - strainCnt[i]);
            if (adiff > eps) {
                converged = false;
            }
            maxDiff = (adiff > maxDiff) ? adiff : maxDiff;
            strainCnt[i] = newStrainCnt[i];
            newStrainCnt[i] = 0.0;
        }

        if (std::abs(readCntValidator - readCnt) > 10) {
            logger->error("original: {}, current : {}, diff : {}", readCnt,
                          readCntValidator, std::abs(readCntValidator - readCnt));
        }
        if (cntr > 0 and cntr % 100 == 0) {
            logger->info("max diff : {}, readCnt : {}", maxDiff, static_cast<uint64_t>(readCntValidator));
        }
    }
    logger->info("iterator cnt: {}", cntr);
    /*auto &eqvec1 = eqb.eqVec();
    std::vector<bool> eqRemained;
    eqRemained.resize(eqvec1.size());
    for (auto kk=0; kk<eqRemained.size(); kk++)
        eqRemained[kk] = false;
    uint32_t kk{0};
    for (auto & eqc: eqvec1) {
        auto &tg = eqc.first;
        auto &v = eqc.second;
        auto csize = v.weights.size();
        for (size_t readMappingCntr = 0; readMappingCntr < csize; ++readMappingCntr) {
            auto &tgt = tg.tgts[readMappingCntr];
            if (strainValid[tgt]) {
                eqRemained[kk] = true;
            }
        }
        kk++;
    }
    kk=0;
    for (auto tf: eqRemained) {
        if (!tf) {
            std::cerr << "so bad eq is removed: " << kk << "\n";
        }
        kk++;
    }*/
    // We have done the EM in the space of sequence / reference IDs
    // but we need to output results in terms of taxa IDs.  Here, we 
    // will map our reference IDs back to taxa IDs, and put the resulting
    // computed abundances in the "strain" member variable.
    decltype(strain) outputMap;
    //spp::sparse_hash_map<uint64_t, double> outputMap;
    outputMap.reserve(strain.size());
    double finalReadCnt{0}, numOfValids{0};
    std::stringstream ss;
    ss << "\n\ndone with EM\n\n";
    std::cerr << ss.str();
    for (auto &kv : strain) {
        finalReadCnt += strainCnt[kv.first];
        numOfValids += strainValid[kv.first];
        util::incLoop(outputMap[seqToTaxMap[kv.first]], strainValid[kv.first] ? strainCnt[kv.first] : 0);
        if (fileReader.refName(kv.first) == "ENST00000495701.1|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000099595.1|CEP104-211|CEP104|485|retained_intron|" or
            fileReader.refName(kv.first) == "ENST00000460038.5|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000009813.1|CEP104-206|CEP104|696|processed_transcript|" or
                fileReader.refName(kv.first) == "ENST00000494653.5|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000009749.1|CEP104-209|CEP104|2174|retained_intron|" or
                fileReader.refName(kv.first) == "ENST00000378230.7|ENSG00000116198.12|OTTHUMG00000003507.5|OTTHUMT00000009747.3|CEP104-202|CEP104|6424|protein_coding|") {
            std::stringstream ss;
            ss << kv.first << " " << fileReader.refName(kv.first) << " " << strainCnt[kv.first] << " " << seqToTaxMap[kv.first] << " " << outputMap[seqToTaxMap[kv.first]]<< "\n";
            std::cerr << ss.str();
        }
    }
    logger->info("Final Reference-level read cnt: {}, # of valid refs: {}",
                 static_cast<uint64_t >(finalReadCnt),
                 static_cast<uint64_t>(numOfValids));
    // Until here, strain map was actually holding refids as key, but after swap it'll be holding strain taxids
    std::swap(strain, outputMap);

    return cntr < maxIter;
}

template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::serialize(std::string &output_filename) {
    logger->info("Write results into the file: {}", output_filename);
    logger->info("# of strains: {}", strain.size());
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\n";
    spp::sparse_hash_map<uint64_t, double> validTaxa;
    for (auto &kv : strain) {
        if (taxaNodeMap.find(kv.first) != taxaNodeMap.end()) {
            TaxaNode *walker = &taxaNodeMap[kv.first];
            while (!walker->isRoot() && walker->getRank() != pruningLevel) {
                walker = &taxaNodeMap[walker->getParentId()];
                if (walker->getId() > 18000000000000) std::exit(1);
            }
            if (!walker->isRoot()) {
                if (validTaxa.find(walker->getId()) == validTaxa.end()) {
                    validTaxa[walker->getId()] = kv.second;
                } else {
                    validTaxa[walker->getId()] += kv.second;
                }
            }
        } else {
            std::cerr << "taxa not found: " << kv.first << "\n";
        }
    }
    double finalReadCnt{0.0};
    for (auto &kv : validTaxa) {
        if (kv.second != 0)
            ofile << kv.first << "\t"
                  << TaxaNode::rank2str(taxaNodeMap[kv.first].getRank())
                  << "\t" << kv.second << "\n";
        finalReadCnt += kv.second;
    }
    ofile.close();
    logger->info("Final reported read count by Cedar: {}", static_cast<uint64_t >(finalReadCnt));
    std::ofstream covOfile(output_filename + ".coverage");
    for (auto &kv: cov) {
        covOfile << kv.first << "\t" << kv.second << "\n";
    }
    covOfile.close();

}

template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::serializeFlat(std::string &output_filename) {
    logger->info("[FlatAbund]");
    // validate final count:
    uint64_t finalMappedReadCnt = 0;
    for (auto &kv : strain) {
        finalMappedReadCnt += kv.second;
    }
    logger->info("Before writing results in the file, total # of mapped reads is {}", finalMappedReadCnt);
    logger->info("Write results in the file: {}", output_filename);
    std::ofstream ofile(output_filename);
    ofile << "taxaId\ttaxaRank\tcount\tcoverage\n";
//    std::cerr << "NUMREFS: " << fileReader.numRefs() << "\n";
    for (uint32_t i = 0; i < fileReader.numRefs(); ++i) {
        //for (auto& kv : strain) {
        auto it = strain.find(i);
        double abund = 0.0;
        if (it != strain.end()) {
            abund = it->second;
        }
        ofile << fileReader.refName(i) << "\t"
              << "flat"
              << "\t" << abund << "\t" << strain_coverage[i] << "\n";
//        std::cerr << fileReader.refName(i) << " " << strain_coverage[i] << "\n";
    }
    ofile.close();
}

template<class ReaderType, class FileReaderType>
void Cedar<ReaderType, FileReaderType>::run(std::string mapperOutput_filename,
                            bool requireConcordance,
                            size_t maxIter,
                            double eps,
                            double minCnt,
                            std::string &output_filename,
                            bool onlyUniq,
                            bool onlyPerf,
                            uint32_t segmentSize,
                            uint32_t rangeFactorizationBins,
                            uint32_t numThreads) {
    loadMappingInfo(mapperOutput_filename, requireConcordance, onlyUniq, onlyPerf,
                    segmentSize, rangeFactorizationBins, numThreads);
    bool verbose = true;
    basicEM(maxIter, eps, minCnt, numThreads, verbose);
    logger->info("serialize to ", output_filename);
    if (!flatAbund) {
        serialize(output_filename);
    } else {
        serializeFlat(output_filename);
    }
    //std::cout << "I guess that's it\n";
}

template
class Cedar<PuffMappingReader, PAMReader>;

template
class Cedar<SAMReader, SAMFileReader>;

/**
 * "How to run" example:
 * make Pufferfish!
 * In the Pufferfish build directory run the following command:
 * /usr/bin/time src/cedar 
 * -t /mnt/scratch2/avi/meta-map/kraken/KrakenDB/taxonomy/nodes.dmp  
 * -s /mnt/scratch2/avi/meta-map/kraken/KrakenDB/seqid2taxid.map 
 * -m /mnt/scratch2/avi/meta-map/kraken/puff/dmps/HC1.dmp 
 * -o HC1.out 
 * -l genus (optional)
 * -f 0.8 (optional)
 **/
int main(int argc, char *argv[]) {
    (void) argc;
    using namespace clipp;
    CedarOpts kopts;
    bool showHelp{false};

    auto checkLevel = [](const char *lin) -> void {
        std::string l(lin);
        std::unordered_set<std::string> valid{"species", "genus", "family", "order", "class", "phylum"};
        if (valid.find(l) == valid.end()) {
            std::string s = "The level " + l + " is not valid.";
            throw std::range_error(s);
        }
    };

    auto cli = (
            (required("--flat").set(kopts.flatAbund, true) % "estimate flat abundance (i.e. there is no taxonomy given)"
             | (
                     required("--taxtree", "-x") &
                     value("taxtree", kopts.taxonomyTree_filename) % "path to the taxonomy tree file",
                             required("--seq2taxa", "-s") &
                             value("seq2taxa", kopts.refId2TaxId_filename) % "path to the refId 2 taxId file "
             )),
                    ((required("--puffMapperOut", "-p").set(kopts.isPuffOut, true) &
                      value("mapout", kopts.mapperOutput_filename) % "path to the pufferfish mapper output file")
                     |
                     (required("--sam").set(kopts.isSAM, true) &
                      value("mapout", kopts.mapperOutput_filename) % "path to the SAM file")
                    ),
                    required("--output", "-o") &
                    value("output", kopts.output_filename) % "path to the output file to write results",
                    option("--maxIter", "-i") &
                    value("iter", kopts.maxIter) % "maximum number of EM iteratons (default : 1000)",
                    option("--eps", "-e") & value("eps", kopts.eps) % "epsilon for EM convergence (default : 0.001)",
                    option("--threads", "-t") & value("threads", kopts.numThreads) % "number of threads to use",
                    option("--rangeFactorizationBins") & value("rangeFactorizationBins", kopts.rangeFactorizationBins) %
                                                         "Number of bins for range factorization (default : 4)",
                    option("--minCnt", "-c") & value("minCnt", kopts.minCnt) %
                                               "minimum count for keeping a reference with count greater than that (default : 0)",
                    option("--level", "-l") & value("level", kopts.level).call(checkLevel) %
                                              "choose between (species, genus, family, order, class, phylum). (default : species)",
                    option("--filter", "-f") & value("filter", kopts.filterThreshold) %
                                               "choose the threshold [0,1] below which to filter out mappings (default : no filter)",
                    option("--noDiscordant").set(kopts.requireConcordance, true) %
                    "ignore orphans for paired end reads",
                    option("--unique").set(kopts.onlyUniq, true) % "report abundance based on unique reads",
                    option("--perfect").set(kopts.onlyPerfect, true) %
                    "report abundance based on perfect reads (unique and with complete coverage)",
                    option("--help", "-h").set(showHelp, true) % "show help",
                    option("-v", "--version").call([] { std::cout << "version 0.1.0\n\n"; }).doc("show version")
    );
    //Multithreaded console logger(with color support)
    auto console = spdlog::stderr_color_mt("console");

    decltype(parse(argc, argv, cli)) res;
    try {
        res = parse(argc, argv, cli);
        if (showHelp) {
            std::cout << make_man_page(cli, "cedar");
            return 0;
        }
    } catch (std::exception &e) {
        std::cout << "\n\nparsing command line failed with exception: " << e.what() << "\n";
        std::cout << "\n\n";
        std::cout << make_man_page(cli, "cedar");
        return 1;
    }

    if (res) {
        if (kopts.isSAM) {
            Cedar<SAMReader, SAMFileReader> cedar(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level,
                                   kopts.filterThreshold, kopts.flatAbund, console);
            cedar.run(kopts.mapperOutput_filename,
                      kopts.requireConcordance,
                      kopts.maxIter,
                      kopts.eps,
                      kopts.minCnt,
                      kopts.output_filename,
                      kopts.onlyUniq,
                      kopts.onlyPerfect,
                      kopts.segmentSize,
                      kopts.rangeFactorizationBins,
                      kopts.numThreads);
        } else {
            Cedar<PuffMappingReader, PAMReader> cedar(kopts.taxonomyTree_filename, kopts.refId2TaxId_filename, kopts.level,
                                           kopts.filterThreshold, kopts.flatAbund, console);
            cedar.run(kopts.mapperOutput_filename,
                      kopts.requireConcordance,
                      kopts.maxIter,
                      kopts.eps,
                      kopts.minCnt,
                      kopts.output_filename,
                      kopts.onlyUniq,
                      kopts.onlyPerfect,
                      kopts.segmentSize,
                      kopts.rangeFactorizationBins,
                      kopts.numThreads);
        }
        return 0;
    } else {
        std::cout << usage_lines(cli, "cedar") << '\n';
        return 1;
    }
}

