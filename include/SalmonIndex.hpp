#ifndef __SALMON_INDEX_HPP__
#define __SALMON_INDEX_HPP__

extern "C" {
#include "bwa.h"
#include "bwamem.h"
#include "kvec.h"
#include "utils.h"
}

#include <memory>

#include <boost/filesystem.hpp>
#include <boost/range/irange.hpp>

#include "format.h"
#include "spdlog/spdlog.h"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"

#include "BWAUtils.hpp"
#include "SalmonConfig.hpp"
#include "SalmonIndexVersionInfo.hpp"
#include "KmerIntervalMap.hpp"

extern "C" {
int bwa_index(int argc, char* argv[]);
}

class SalmonIndex{
        public:
            SalmonIndex(std::shared_ptr<spdlog::logger>& logger) : 
                loaded_(false), versionInfo_(0, false, 0), logger_(logger) {}

            ~SalmonIndex() {
                if (idx_) { bwa_idx_destroy(idx_); }
            }
           
            void load(const boost::filesystem::path& indexDir) { 
                namespace bfs = boost::filesystem;

                // Check if version file exists and, if so, read it.
                boost::filesystem::path versionPath = indexDir / "versionInfo.json";
                versionInfo_.load(versionPath); 
                if (versionInfo_.indexVersion() == 0) {
                    fmt::MemoryWriter infostr;
                    infostr << "Error: The index version file " << versionPath.string() 
                            << " doesn't seem to exist.  Please try re-building the salmon "
                               "index.";
                    throw std::invalid_argument(infostr.str());
                }
                // Check index version compatibility here 

                if (versionInfo_.hasAuxKmerIndex()) {
                    // Read the aux index
                    logger_->info("Loading auxiliary index");
                    bfs::path auxIdxFile = indexDir / "aux.idx";
                    auxIdx_.setK(versionInfo_.auxKmerLength());
                    auxIdx_.load(auxIdxFile);
                    logger_->info("Auxiliary index contained {} k-mers", auxIdx_.size());
                    logger_->info("done");
                }

                logger_->info("Loading BWA index");
                // Read the actual BWA index
                { // mem-based
                    boost::filesystem::path indexPath = indexDir / "bwaidx";
                    //if ((idx_ = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_BWT|BWA_IDX_BNS|BWA_IDX_PAC)) == 0) {
                    if ((idx_ = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_ALL)) == 0) {
                        fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                        fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                        std::exit(1);
                    }
                }
                logger_->info("done");

                loaded_ = true;
            }
            
            bool buildAux_(boost::filesystem::path indexDir, uint32_t k) {
                       namespace bfs = boost::filesystem;
                       
                       bfs::path indexPath = indexDir / "bwaidx";
                       // Load the bwa index
                       {
                           logger_->info("Reading BWT index from file");
                           if ((idx_ = bwa_idx_load(indexPath.string().c_str(), BWA_IDX_BWT|BWA_IDX_BNS|BWA_IDX_PAC)) == 0) {
                               logger_->error("Couldn't open index [{}] --- ", indexPath);
                               logger_->error("Please make sure that 'salmon index' has been run successfully");
                               std::exit(1);
                           }
                       }

                       auxIdx_.setK(k);
                       size_t numRecords = idx_->bns->n_seqs;
                       { // Load transcripts from file
                          logger_->info("Index contained {} targets; streaming through them", numRecords);
                          for (auto i : boost::irange(size_t(0), numRecords)) {
                              char* name = idx_->bns->anns[i].name;
                              uint32_t len = idx_->bns->anns[i].len;
                              uint8_t* rseq = nullptr;
                              int64_t tstart, tend, compLen, l_pac = idx_->bns->l_pac;
                              tstart  = idx_->bns->anns[i].offset;
                              tend = tstart + len;
                              rseq = bns_get_seq(l_pac, idx_->pac, tstart, tend, &compLen);
                              if (compLen != len) {
                                  fmt::print(stderr,
                                          "For transcript {}, stored length ({}) != computed length ({}) --- index may be corrupt. exiting\n",
                                          name, compLen, len);
                                  std::exit(1);
                              }
                              if (len < k) { continue; }
                              for (int s = 0; s < len - k + 1; ++s) {
                                  bwtintv_t resInterval;
                                  KmerKey key(&(rseq[s]), k);
                                  if (!auxIdx_.hasKmer(key)) {
                                      bool foundInterval = bwautils::getIntervalForKmer(idx_->bwt, k, &(rseq[s]), resInterval);
                                      // If we found the interval for this k-mer
                                      if (foundInterval) {
                                          // If the interval for this k-mer isn't already in the hash
                                          // then put it in the hash
                                          auxIdx_[key] = resInterval;
                                      }
                                  }
                              }
                          }
                          // Since we have the de-coded reference sequences, we no longer need
                          // the encoded sequences, so free them.
                          free(idx_->pac); idx_->pac = nullptr;
                          // ====== Done streaming through transcripts 
                       }
                      
                       bfs::path auxIndexFile = indexDir / "aux.idx"; 
                       auxIdx_.save(auxIndexFile);
                       return true;
            }


            bool build(boost::filesystem::path indexDir, 
                       std::vector<const char*>& bwaArgVec,
                       uint32_t k) {
                namespace bfs = boost::filesystem;
                char* bwaArgv[] = { const_cast<char*>(bwaArgVec[0]),
                    const_cast<char*>(bwaArgVec[1]),
                    const_cast<char*>(bwaArgVec[2]),
                    const_cast<char*>(bwaArgVec[3]),
                    const_cast<char*>(bwaArgVec[4]),
                    const_cast<char*>(bwaArgVec[5]) };
                int bwaArgc = 6;
                int ret = bwa_index(bwaArgc, bwaArgv);

                bool buildAux = (k > 0);
                if (buildAux) {
                    buildAux_(indexDir, k);
                }

                bfs::path versionFile = indexDir / "versionInfo.json";
                versionInfo_.indexVersion(salmon::indexVersion); 
                versionInfo_.hasAuxKmerIndex(buildAux);
                versionInfo_.auxKmerLength(k);
                versionInfo_.save(versionFile);
                return (ret == 0);
            }

            bool loaded() { return loaded_; }
            bwaidx_t* bwaIndex() { return idx_; }

            bool hasAuxKmerIndex() { return versionInfo_.hasAuxKmerIndex(); }
            KmerIntervalMap& auxIndex() { return auxIdx_; }
        private:
          bool loaded_;
          SalmonIndexVersionInfo versionInfo_;
          bwaidx_t *idx_{nullptr};
          KmerIntervalMap auxIdx_;
          std::shared_ptr<spdlog::logger> logger_;
};

#endif //__SALMON_INDEX_HPP
