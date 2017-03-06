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

#include "spdlog/spdlog.h"
#include "cereal/archives/json.hpp"
#include "cereal/types/vector.hpp"

#include "BooMap.hpp"
#include "FrugalBooMap.hpp"
#include "RapMapSAIndex.hpp"
#include "IndexHeader.hpp"
#include "BWAUtils.hpp"
#include "SalmonConfig.hpp"
#include "SalmonIndexVersionInfo.hpp"
#include "KmerIntervalMap.hpp"

extern "C" {
int bwa_index(int argc, char* argv[]);
}
// declaration of quasi index function
int rapMapSAIndex(int argc, char* argv[]);

template <typename IndexT> 
using DenseHash = spp::sparse_hash_map<uint64_t, 
                                         rapmap::utils::SAInterval<IndexT>, 
                                         rapmap::utils::KmerKeyHasher>;
template <typename IndexT> 
using PerfectHash = FrugalBooMap<uint64_t, rapmap::utils::SAInterval<IndexT>>;

class SalmonIndex{
        public:
            SalmonIndex(std::shared_ptr<spdlog::logger>& logger, SalmonIndexType indexType) :
                loaded_(false), versionInfo_(0, false, 0, indexType), logger_(logger), seqHash_(""), nameHash_("") {}

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

                auto indexType = versionInfo_.indexType();
                // Load the appropriate index type
                if (indexType == SalmonIndexType::FMD) {
                    loadFMDIndex_(indexDir);
                } else {
                    loadQuasiIndex_(indexDir);
                }

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
                                  logger_->error(
                                          "For transcript {}, stored length ({}) != computed length ({}) --- index may be corrupt. exiting\n",
                                          name, compLen, len);
                                  std::exit(1);
                              }
                              if (len < k) { continue; }
                              for (uint32_t s = 0; s < len - k + 1; ++s) {
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
                       std::vector<std::string>& argVec,
                       uint32_t k) {
                namespace bfs = boost::filesystem;
                switch (versionInfo_.indexType()) {
                    case SalmonIndexType::QUASI:
                        return buildQuasiIndex_(indexDir, argVec, k);
                    case SalmonIndexType::FMD:
                        return buildFMDIndex_(indexDir, argVec, k);
                    default:
                        logger_->warn("Unexpected index type; cannot build");
                        return false;
                }
            }

            bool loaded() { return loaded_; }
            bwaidx_t* bwaIndex() { return idx_; }

            bool is64BitQuasi() { return largeQuasi_; }
            bool isPerfectHashQuasi() { return perfectHashQuasi_;} 

            RapMapSAIndex<int32_t, DenseHash<int32_t>>* quasiIndex32() { return quasiIndex32_.get(); }
            RapMapSAIndex<int64_t, DenseHash<int64_t>>* quasiIndex64() { return quasiIndex64_.get(); }

            RapMapSAIndex<int32_t, PerfectHash<int32_t>>* quasiIndexPerfectHash32() { return quasiIndexPerfectHash32_.get(); }
            RapMapSAIndex<int64_t, PerfectHash<int64_t>>* quasiIndexPerfectHash64() { return quasiIndexPerfectHash64_.get(); }

            bool hasAuxKmerIndex() { return versionInfo_.hasAuxKmerIndex(); }
            KmerIntervalMap& auxIndex() { return auxIdx_; }

            SalmonIndexType indexType() { return versionInfo_.indexType(); }

	    const char* transcriptomeSeq() {
	      if (loaded_) {
		if (is64BitQuasi()) {
		  return quasiIndex64_->seq.c_str();
		} else {
		  return quasiIndex32_->seq.c_str();
		}
	      } else {
		return nullptr;
	      }
	    }

	    uint64_t transcriptOffset(uint64_t id) {
	      if (loaded_) {
		if (is64BitQuasi()) {
		  return quasiIndex64_->txpOffsets[id];
		} else {
		  return quasiIndex32_->txpOffsets[id];
		}
	      } else {
		return std::numeric_limits<uint64_t>::max();
	      }
	    }

	std::string seqHash() const { return seqHash_; }
	std::string nameHash() const { return nameHash_; }

        private:
            bool buildFMDIndex_(boost::filesystem::path indexDir,
                                std::vector<std::string>& bwaArgVec,
                                uint32_t k) {
                namespace bfs = boost::filesystem;
                char* bwaArgv[] = {
                    const_cast<char*>(bwaArgVec[0].c_str()),
                    const_cast<char*>(bwaArgVec[1].c_str()),
                    const_cast<char*>(bwaArgVec[2].c_str()),
                    const_cast<char*>(bwaArgVec[3].c_str()),
                    const_cast<char*>(bwaArgVec[4].c_str()),
                    const_cast<char*>(bwaArgVec[5].c_str()) };
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

            bool buildQuasiIndex_(boost::filesystem::path indexDir,
                                  std::vector<std::string>& quasiArgVec,
                                  uint32_t k) {
                namespace bfs = boost::filesystem;
		int quasiArgc = static_cast<int>(quasiArgVec.size());
		char** quasiArgv = new char*[quasiArgc];
		for (size_t i = 0; i < quasiArgc; ++i) {
		  auto& arg = quasiArgVec[i];
		  quasiArgv[i] = new char[arg.size() + 1];
		  std::strcpy(quasiArgv[i], arg.c_str());
		}

                int ret = rapMapSAIndex(quasiArgc, quasiArgv);

                bfs::path versionFile = indexDir / "versionInfo.json";
                versionInfo_.indexVersion(salmon::indexVersion);
                versionInfo_.hasAuxKmerIndex(false);
                versionInfo_.auxKmerLength(k);
                versionInfo_.indexType(SalmonIndexType::QUASI);
                versionInfo_.save(versionFile);
	
		// Free the memory used for the arg vector
		for (size_t i = 0; i < quasiArgc; ++i) {
		  delete quasiArgv[i];
		}
		delete [] quasiArgv;
	
		return (ret == 0);
            }

          bool loadFMDIndex_(const boost::filesystem::path& indexDir) {
              namespace bfs = boost::filesystem;
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
                      logger_->error("Couldn't open index [{}] --- ", indexPath);
                      logger_->error("Please make sure that 'salmon index' has been run successfully");
                      std::exit(1);
                  }
              }
              logger_->info("done");
              return true;
          }

          bool loadQuasiIndex_(const boost::filesystem::path& indexDir) {
              namespace bfs = boost::filesystem;
              logger_->info("Loading Quasi index");
              // Read the actual Quasi index
              { // quasi-based
                  boost::filesystem::path indexPath = indexDir;
                  std::string indexStr = indexDir.string();
                  if (indexStr.back() != '/') { indexStr.push_back('/'); }

                  IndexHeader h;
                  std::ifstream indexStream(indexStr + "header.json");
                  {
                    cereal::JSONInputArchive ar(indexStream);
                    ar(h);
                  }
                  indexStream.close();

                  if (h.version() != salmon::requiredQuasiIndexVersion) {
		    logger_->critical("I found a quasi-index with version {}, but I require {}. "
				      "Please re-index the reference.", h.version(), salmon::requiredQuasiIndexVersion);
                    std::exit(1);
                  }
                  if (h.indexType() != IndexType::QUASI) {
                    logger_->critical("The index {} does not appear to be of the "
                                      "appropriate type (quasi)", indexStr);
                    std::exit(1);
                  }

		  seqHash_ = h.seqHash();
		  nameHash_ = h.nameHash();

                  // Is the quasi-index using a perfect hash
                  perfectHashQuasi_ = h.perfectHash();

                  if (h.bigSA()) {
                    largeQuasi_ = true;
                    logger_->info("Loading 64-bit quasi index");
                    if (perfectHashQuasi_) {
                        quasiIndexPerfectHash64_.reset(new RapMapSAIndex<int64_t, PerfectHash<int64_t>>);
                        if (!quasiIndexPerfectHash64_->load(indexStr)) {
                            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                            fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                            std::exit(1);
                        }
                    } else {
                        quasiIndex64_.reset(new RapMapSAIndex<int64_t, DenseHash<int64_t>>);
                        if (!quasiIndex64_->load(indexStr)) {
                            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                            fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                            std::exit(1);
                        }
                    }
                  } else { // 32-bit index
                    logger_->info("Loading 32-bit quasi index");
                    
                    if (perfectHashQuasi_) {
                        quasiIndexPerfectHash32_.reset(new RapMapSAIndex<int32_t, PerfectHash<int32_t>>);
                        if (!quasiIndexPerfectHash32_->load(indexStr)) {
                            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                            fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                            std::exit(1);
                        }
                    } else {
                        quasiIndex32_.reset(new RapMapSAIndex<int32_t, DenseHash<int32_t>>);
                        if (!quasiIndex32_->load(indexStr)) {
                            fmt::print(stderr, "Couldn't open index [{}] --- ", indexPath);
                            fmt::print(stderr, "Please make sure that 'salmon index' has been run successfully\n");
                            std::exit(1);
                        }
                    }
                  }
              }
              logger_->info("done");
              return true;
          }


          bool loaded_;
          SalmonIndexVersionInfo versionInfo_;
          // Can't think of a generally better way to do this now
          // without making the entire code-base look crazy
          bool largeQuasi_{false};
          bool perfectHashQuasi_{false};

          std::unique_ptr<RapMapSAIndex<int32_t, DenseHash<int32_t>>> quasiIndex32_{nullptr};
          std::unique_ptr<RapMapSAIndex<int64_t, DenseHash<int64_t>>> quasiIndex64_{nullptr};

          std::unique_ptr<RapMapSAIndex<int32_t, PerfectHash<int32_t>>> quasiIndexPerfectHash32_{nullptr};
          std::unique_ptr<RapMapSAIndex<int64_t, PerfectHash<int64_t>>> quasiIndexPerfectHash64_{nullptr};

          bwaidx_t *idx_{nullptr};
          KmerIntervalMap auxIdx_;
          std::shared_ptr<spdlog::logger> logger_;
	  std::string seqHash_;
	  std::string nameHash_;
};

#endif //__SALMON_INDEX_HPP
